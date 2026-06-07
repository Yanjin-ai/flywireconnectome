# Deploying the interactive explorer

The Streamlit explorer (`src/explorer_app.py`) runs entirely from committed
artifacts (`network_enriched.csv`, `results/*.json`) — **no 300 MB data download
is needed**, so it deploys to Streamlit Community Cloud for free.

## One-click deploy (Streamlit Community Cloud)

1. Click the badge (or this link), which pre-fills the deploy form:

   [![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/deploy?repository=Yanjin-ai/flywireconnectome&branch=main&mainModule=src/explorer_app.py)

   <https://share.streamlit.io/deploy?repository=Yanjin-ai/flywireconnectome&branch=main&mainModule=src/explorer_app.py>

2. Sign in with the GitHub account that owns the repo (one-time OAuth — this
   step can only be done by you; it can't be automated).
3. Confirm: **Repository** `Yanjin-ai/flywireconnectome`, **Branch** `main`,
   **Main file path** `src/explorer_app.py`. **Leave the Python version at the
   default** that Streamlit offers (do *not* manually switch it).
   > Known quirk: changing the Python dropdown (e.g. to 3.12) makes the deploy
   > form re-validate into a bad state and falsely show "This repository/branch/
   > file does not exist". At the default version the fields validate correctly.
   > `requirements.txt` now uses version *floors*, so any offered Python works.
4. Click **Deploy**. First build installs `requirements.txt` (~1–2 min).

Your app will be live at a URL like
`https://flywireconnectome-<hash>.streamlit.app` (you can set a custom subdomain
in the app settings, e.g. `https://flywire-mcis.streamlit.app`).

## Run locally

```bash
pip install -e ".[test]"        # or: pip install -r requirements.txt
streamlit run src/explorer_app.py
```

## What the explorer shows

- Headline metrics (N, conserved edges, 7.4× edge conservation vs degree-null,
  descending enrichment).
- **Conserved subgraph** tab — the conserved directed edges by cell type.
- **Conservation track** tab — per-neuron conservation z-scores (top neurons +
  the ranked track curve).
- **Composition & download** tab — class / neurotransmitter breakdown and a
  filtered-circuit CSV download.

The app reads only files already in the repo, so the deployed demo always
matches the committed results.
