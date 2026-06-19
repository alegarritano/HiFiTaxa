#!/usr/bin/env python3
"""Append GitHub clone + view traffic to CSVs, deduped by date.

GitHub only keeps 14 days of traffic. Running this on any interval shorter than
14 days turns that rolling window into permanent, gap-free daily history. Re-runs
are idempotent: a date already recorded is just overwritten with the latest value.

Usage:
    TRAFFIC_TOKEN=<token> python3 bin/log_traffic.py <owner/repo> [out_dir]

The token needs permission to read this repo's traffic:
  - fine-grained PAT with "Administration: Read" on the repo, or
  - classic PAT with the "repo" scope.
(GITHUB_TOKEN also works as the env var name; the repo's default Actions token
 does NOT have traffic access, so a PAT is required.)

Writes <out_dir>/clones.csv and <out_dir>/views.csv (default out_dir: traffic).
"""
import csv, json, os, sys, urllib.request


def fetch(repo, kind, token):
    url = f"https://api.github.com/repos/{repo}/traffic/{kind}"
    req = urllib.request.Request(url, headers={
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    })
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def to_rows(payload, list_key):
    return [
        {"date": e["timestamp"][:10], "count": e["count"], "uniques": e["uniques"]}
        for e in payload.get(list_key, [])
    ]


def upsert_csv(path, rows):
    """Merge rows (keyed by date) into an existing CSV, keeping full history."""
    by_date = {}
    if os.path.exists(path):
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                by_date[row["date"]] = row
    for row in rows:
        by_date[row["date"]] = row
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["date", "count", "uniques"])
        w.writeheader()
        for d in sorted(by_date):
            w.writerow(by_date[d])
    return len(by_date)


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: log_traffic.py <owner/repo> [out_dir]")
    repo = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "traffic"
    token = os.environ.get("TRAFFIC_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if not token:
        sys.exit("set TRAFFIC_TOKEN to a PAT that can read this repo's traffic")

    clones = fetch(repo, "clones", token)
    views = fetch(repo, "views", token)
    n_c = upsert_csv(os.path.join(out_dir, "clones.csv"), to_rows(clones, "clones"))
    n_v = upsert_csv(os.path.join(out_dir, "views.csv"), to_rows(views, "views"))
    print(f"[traffic] clones.csv: {n_c} day(s) total; views.csv: {n_v} day(s) total")
    print(f"[traffic] latest 14-day window: "
          f"{clones.get('count', 0)} clones ({clones.get('uniques', 0)} unique), "
          f"{views.get('count', 0)} views ({views.get('uniques', 0)} unique)")


if __name__ == "__main__":
    main()
