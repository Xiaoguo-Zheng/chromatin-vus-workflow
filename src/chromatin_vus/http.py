"""Small auditable HTTP client with bounded retries and content-addressed caches."""
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from datetime import datetime, timezone
import hashlib
import json
import time


class ProviderError(RuntimeError):
    pass


def utc_now():
    return datetime.now(timezone.utc).isoformat()


class Client:
    def __init__(self, cache, offline=False, refresh=False, timeout=30, ttl_days=7):
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.offline, self.refresh, self.timeout, self.ttl_days = offline, refresh, timeout, ttl_days
        self.events = []
        self.last_request = 0.0

    def get(self, url):
        key = hashlib.sha256(url.encode()).hexdigest()
        path = self.cache / (key + ".json")
        if path.exists() and not self.refresh:
            try:
                stored = json.loads(path.read_text(encoding="utf-8"))
                serialized = json.dumps(stored["payload"], sort_keys=True, ensure_ascii=False).encode()
                if stored["url"] != url or hashlib.sha256(serialized).hexdigest() != stored["sha256"]:
                    raise ValueError("cache integrity failure")
                age = time.time() - datetime.fromisoformat(stored["retrieved_at"]).timestamp()
                if self.offline or age <= self.ttl_days * 86400:
                    self.events.append({k: stored[k] for k in ["url", "retrieved_at", "sha256"]} | {"cache": True, "stale": age > self.ttl_days * 86400})
                    return stored["payload"]
            except (ValueError, KeyError, TypeError) as exc:
                if self.offline:
                    raise ProviderError(f"Invalid cached response: {exc}") from exc
        if self.offline:
            raise ProviderError("No valid cached response is available in offline mode.")
        error = None
        for attempt in range(3):
            time.sleep(max(0, 0.2 - (time.monotonic() - self.last_request)))
            self.last_request = time.monotonic()
            try:
                with urlopen(Request(url, headers={"Accept": "application/json", "User-Agent": "chromatin-vus-workflow/0.1.0"}), timeout=self.timeout) as response:
                    data = response.read(16_000_001)
                if len(data) > 16_000_000:
                    raise ProviderError("Provider response exceeds 16 MB.")
                obj = json.loads(data)
                serialized = json.dumps(obj, sort_keys=True, ensure_ascii=False).encode()
                stored = {"url": url, "retrieved_at": utc_now(), "sha256": hashlib.sha256(serialized).hexdigest(), "payload": obj}
                temporary = path.with_suffix(".tmp")
                temporary.write_text(json.dumps(stored, ensure_ascii=False), encoding="utf-8")
                temporary.replace(path)
                self.events.append({k: stored[k] for k in ["url", "retrieved_at", "sha256"]} | {"cache": False, "stale": False})
                return obj
            except HTTPError as exc:
                error = f"HTTP {exc.code}"
                if exc.code == 404:
                    self.events.append({"url": url, "retrieved_at": utc_now(), "error": error})
                    raise ProviderError(error) from exc
                if exc.code not in {429, 500, 502, 503, 504}:
                    break
                try:
                    delay = float(exc.headers.get("Retry-After", 2 ** attempt))
                except (TypeError, ValueError):
                    delay = 2 ** attempt
                if delay > 10:
                    error += "; provider requests a longer retry delay; rerun later"
                    break
                if attempt < 2:
                    time.sleep(max(0, delay))
            except (URLError, TimeoutError, ValueError, OSError) as exc:
                error = f"{type(exc).__name__}: {str(exc)[:180]}"
                if attempt < 2:
                    time.sleep(2 ** attempt)
        self.events.append({"url": url, "retrieved_at": utc_now(), "error": error})
        raise ProviderError(error or "Provider unavailable")
