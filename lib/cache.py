"""
cache.py
Networking with retry/backoff plus a TTL disk cache.

v2.x issued one bare ``requests.get`` per call with no retry, no cache and no
offline mode.  On a firewalled host every single call failed; on a working host
``cpd_id_map`` issued one HTTP round-trip *per compound*.

Public API
----------
  cache_dir      : resolve (and create) the cache directory
  http_get       : cached GET returning text
  http_get_bytes : cached GET returning bytes
  http_post_json : POST returning parsed JSON (not cached)
  set_offline    : force offline mode (cache-only); useful in CI
  is_offline     : query offline mode
  clear_cache    : delete cached payloads
"""

from __future__ import annotations

import hashlib
import json
import os
import time
import warnings
from pathlib import Path
from typing import Any

from .constants import CACHE_TTL_SECONDS, HTTP_BACKOFF, HTTP_RETRIES, HTTP_TIMEOUT
from .errors import NetworkError

_OFFLINE = os.environ.get("PATHVIEW_OFFLINE", "").lower() in ("1", "true", "yes")
_SESSION = None


# ---------------------------------------------------------------------------
# Offline switch
# ---------------------------------------------------------------------------

def set_offline(offline: bool = True) -> None:
    """Force cache-only operation.  No network call will be attempted."""
    global _OFFLINE
    _OFFLINE = bool(offline)


def is_offline() -> bool:
    """True when pathview-plus is running in cache-only mode."""
    return _OFFLINE


# ---------------------------------------------------------------------------
# Cache location
# ---------------------------------------------------------------------------

def cache_dir() -> Path:
    """
    Return the cache directory, creating it if needed.

    Honours ``PATHVIEW_CACHE`` then ``XDG_CACHE_HOME``, else ``~/.cache``.
    """
    env = os.environ.get("PATHVIEW_CACHE")
    if env:
        root = Path(env)
    else:
        xdg = os.environ.get("XDG_CACHE_HOME")
        root = Path(xdg) / "pathview-plus" if xdg else Path.home() / ".cache" / "pathview-plus"
    root.mkdir(parents=True, exist_ok=True)
    return root


def _cache_path(url: str, suffix: str) -> Path:
    key = hashlib.sha256(url.encode()).hexdigest()[:32]
    return cache_dir() / f"{key}{suffix}"


def _read_cache(path: Path, ttl: int) -> bytes | None:
    if not path.exists():
        return None
    if ttl >= 0 and (time.time() - path.stat().st_mtime) > ttl:
        return None
    try:
        return path.read_bytes()
    except OSError:
        return None


def clear_cache() -> int:
    """Delete every cached payload.  Returns the number of files removed."""
    n = 0
    for p in cache_dir().glob("*"):
        if p.is_file():
            p.unlink()
            n += 1
    return n


# ---------------------------------------------------------------------------
# Session with retries
# ---------------------------------------------------------------------------

def _session():
    global _SESSION
    if _SESSION is not None:
        return _SESSION
    import requests
    from requests.adapters import HTTPAdapter

    sess = requests.Session()
    sess.headers.update({
        "User-Agent": "pathview-plus/3.0 (https://github.com/raw-lab/pathview-plus)"
    })
    try:
        from urllib3.util.retry import Retry
        retry = Retry(
            total=HTTP_RETRIES,
            backoff_factor=HTTP_BACKOFF,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset(["GET", "POST"]),
        )
        adapter = HTTPAdapter(max_retries=retry, pool_maxsize=16)
    except Exception:                                    # pragma: no cover
        adapter = HTTPAdapter(max_retries=HTTP_RETRIES, pool_maxsize=16)
    sess.mount("https://", adapter)
    sess.mount("http://", adapter)
    _SESSION = sess
    return sess


# ---------------------------------------------------------------------------
# Cached GET
# ---------------------------------------------------------------------------

def http_get_bytes(
    url: str,
    ttl: int = CACHE_TTL_SECONDS,
    timeout: int = HTTP_TIMEOUT,
    suffix: str = ".bin",
    use_cache: bool = True,
) -> bytes:
    """
    GET *url* and return the body as bytes, using the disk cache when fresh.

    Raises NetworkError when the resource is unavailable and no cached copy
    exists.  A stale cached copy is preferred over a hard failure.
    """
    path = _cache_path(url, suffix)

    if use_cache:
        hit = _read_cache(path, ttl)
        if hit is not None:
            return hit

    if _OFFLINE:
        stale = _read_cache(path, ttl=-1)
        if stale is not None:
            return stale
        raise NetworkError(f"Offline mode and no cached copy of {url}")

    try:
        resp = _session().get(url, timeout=timeout)
        resp.raise_for_status()
        body = resp.content
    except Exception as exc:
        stale = _read_cache(path, ttl=-1)
        if stale is not None:
            warnings.warn(
                f"{url} unreachable ({exc}); using stale cached copy.",
                stacklevel=2,
            )
            return stale
        raise NetworkError(f"Could not fetch {url}: {exc}") from exc

    if use_cache:
        try:
            path.write_bytes(body)
        except OSError:                                   # pragma: no cover
            pass
    return body


def http_get(url: str, **kw: Any) -> str:
    """GET *url* and return the body decoded as UTF-8."""
    kw.setdefault("suffix", ".txt")
    return http_get_bytes(url, **kw).decode("utf-8", errors="replace")


def http_post_json(
    url: str,
    payload: dict,
    timeout: int = HTTP_TIMEOUT,
    ttl: int = CACHE_TTL_SECONDS,
) -> Any:
    """
    POST *payload* as form data and return parsed JSON.

    POST responses are cached on a hash of (url, payload) because the identifier
    mapping services are idempotent and are the dominant cost of a mapping run.
    """
    key = url + "|" + json.dumps(payload, sort_keys=True)
    path = _cache_path(key, ".json")

    hit = _read_cache(path, ttl)
    if hit is not None:
        try:
            return json.loads(hit)
        except json.JSONDecodeError:
            pass

    if _OFFLINE:
        raise NetworkError(f"Offline mode and no cached copy of POST {url}")

    try:
        resp = _session().post(url, data=payload, timeout=timeout)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        raise NetworkError(f"POST {url} failed: {exc}") from exc

    try:
        path.write_bytes(json.dumps(data).encode())
    except OSError:                                       # pragma: no cover
        pass
    return data
