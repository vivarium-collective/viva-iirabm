#!/usr/bin/env bash
# Preview the published read-only workbench LOCALLY.
#
# The bundle built by publish_dashboard.sh uses absolute asset URLs rooted at
# its production base path (/viva-iirabm/workbench/, for GitHub Pages). Opening
# index.html as a file:// therefore 404s every stylesheet/script and renders an
# unstyled page. This script serves the bundle over HTTP at the matching path so
# it renders exactly as it will in production.
#
# Usage: scripts/preview_dashboard.sh [port]
set -euo pipefail
WS_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUNDLE="$WS_ROOT/reports/published/workbench"
PORT="${1:-8913}"
BASE="/viva-iirabm/workbench"

if [ ! -f "$BUNDLE/index.html" ]; then
  echo "No bundle yet. Build it first:"
  echo "  bash scripts/publish_dashboard.sh reports/published/workbench \"$BASE\""
  exit 1
fi

SERVE="$(mktemp -d)"
mkdir -p "$SERVE/viva-iirabm"
ln -s "$BUNDLE" "$SERVE/viva-iirabm/workbench"
URL="http://127.0.0.1:${PORT}${BASE}/"
echo "Serving $BUNDLE at $URL  (Ctrl-C to stop)"
command -v open >/dev/null && ( sleep 1; open "$URL" ) &
cd "$SERVE"
exec python3 -m http.server "$PORT" --bind 127.0.0.1
