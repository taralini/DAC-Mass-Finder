# Antibody–Peptide Mass Matcher (Web)

A simple Streamlit web app to match observed masses (neutral or m/z) to predicted antibody–peptide conjugation stoichiometries.

## Quick Deploy Options

### Option A — Streamlit Community Cloud (easiest)
1. Create a new GitHub repo and add these files.
2. Go to https://share.streamlit.io/ and click **Deploy an app**.
3. Select your repo and set **Main file path** to `app.py`.
4. (Optional) Set **Advanced settings → Python version** to 3.11+.
5. Click **Deploy**.

### Option B — Render (Docker)
1. Push this folder to a GitHub repo.
2. Create a new **Web Service** on https://render.com
3. Choose **Docker** and point to your repo; Render will pick up `render.yaml` or you can specify the Dockerfile directly.
4. Deploy.

### Option C — Fly.io (Docker)
```bash
flyctl launch --no-deploy
# edit fly.toml if needed
flyctl deploy
```

### Option D — Heroku
```bash
heroku create adc-mass-matcher
git push heroku main
heroku ps:scale web=1
heroku open
```

### Option E — Any Docker Host
```bash
docker build -t adc-mass-matcher .
docker run -p 8501:8501 adc-mass-matcher
# open http://localhost:8501
```

## Local Dev
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

## CSV Formats

**peptides.csv**
```
name,peptide_mass,linker_mass,max_copies
PepA,1998.987,1200.123,8
PepB,2999.020,0,4
```

**observed.csv**
```
observed_mass
155800.1
158900.0
```

## Notes
- Set `Δ per conjugation` to account for chemistry-dependent mass changes (e.g., -18.01056).
- Use **ppm** tolerance for high-resolution MS; otherwise **Da** tolerance.
- If your observed data are m/z, enter charge states so the app computes predicted m/z accordingly.
