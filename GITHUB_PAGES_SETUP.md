# Publish Aircraft Design Lab on your GitHub account

This package is a complete static PyScript website. Students only need a current
browser and internet access. They do not need Python, Anaconda or a login.

## Publish using the GitHub website

1. Download `aircraft_design_pyscript.zip` and extract it on your computer.
2. Sign into **your** GitHub account at https://github.com .
3. Create a new **public** repository named `Aircraft-Design`. Creating a separate
   repository leaves your `Helicopter-Design` site unchanged. Initialize the new
   repository with a README so that the `main` branch exists.
4. In the new repository, select **Add file → Upload files**.
5. Open the extracted package folder. Drag its **contents** into the GitHub
   upload area, including the `model`, `examples` and `reference` folders. Keep
   their folder structure. Do not upload the ZIP itself or add another outer
   `aircraft_design_pyscript` folder.
6. Commit the upload. Check that `index.html` is visible at the top level of the
   repository, alongside `app.js`, `style.css`, `browser_runner.py`, `pyscript.json`,
   `desktop_model.zip` and the folders. The package also includes `.nojekyll`.
7. Open **Settings → Pages**. Under **Build and deployment**, set **Source** to
   **Deploy from a branch**, choose **main** and **/(root)**, then **Save**.
8. Wait for GitHub to publish. Settings → Pages will show the actual live URL.
   Publication can take up to 10 minutes. Use the URL GitHub gives you.

If the username is `aeroguruhk` and the repository is exactly `Aircraft-Design`,
the expected address is:

`https://aeroguruhk.github.io/Aircraft-Design/`

That is an expected address, not a claim that it has already been published.
The separate public link provided in our chat is already hosted elsewhere;
GitHub publication requires the steps above in your account.

## How students use the site

1. Choose Fighter jet, Transport jet, GA turboprop or Transport turboprop.
2. Adjust mission, wing or propulsion inputs, or edit/load a complete JSON.
3. Press **Run sizing analysis**.
4. First use downloads Python and numerical packages. Allow a few minutes on
   a slow connection; a desktop or laptop is recommended for the first lesson.
5. View constraint, weight, wing-load and performance tabs.
6. Download the full ZIP of generated plots, CSVs, input JSON and summary.

The initial results are clearly marked **Reference example**. They are previously
verified desktop results, not a new browser calculation. Editing inputs marks
those results stale. A successful fresh run is marked **Browser calculation**.

Calculations run on the student's device. The application does not upload entered
JSON to a calculation server. All computation is client-side; internet is needed
for loading the site, PyScript and scientific packages. It is not an offline bundle.

## Browser support

Designed for current desktop Chrome, Edge, Firefox and Safari with WebAssembly,
JavaScript and workers enabled. No claim is made for older browsers or every
mobile device. School network filters may need to allow `pyscript.net`,
`cdn.jsdelivr.net`, and PyScript's package resources. PyScript is pinned to
2026.7.3, and the Pyodide Python runtime to 0.28.3.

The worker uses `sync_main_only: true`; the Python calculation does not access
the browser DOM. This supports ordinary GitHub Pages hosting without special
cross-origin-isolation headers or a service worker.

## Local preview on Windows / Anaconda

Do not double-click `index.html`; file:// pages cannot fetch the Python modules.
Open Anaconda Prompt in the extracted folder and run:

```text
python -m http.server 8000
```

Then visit `http://localhost:8000/` in your browser. The local Python command only
serves the files; the calculation itself still runs inside the browser. The
published HTTPS page is what you should share with students.

## Editing defaults and updating

- Change example defaults in `examples/*.json`.
- Main-page text is in `index.html`, style in `style.css`, browser interactions
  in `app.js`, runtime configuration in `pyscript.json`.
- Numerical modules in `model/` are the same as the verified desktop V2 model.
- If examples or model code change, regenerate or remove the corresponding
  `reference/` results so reference cards do not imply verification of new inputs.
- Upload replacements and commit; GitHub Pages republishes the selected branch.
- If students see older files after an update, reload the page and clear the
  browser cache if necessary. Do not update runtime versions without retesting.

## What was verified

- All four examples were run in Pyodide 0.28.3's WebAssembly Python engine, with
  NumPy, SciPy and Matplotlib. Gross masses matched desktop values exactly.
- All four runs generated plots and result ZIPs, which were checked as ZIP archives.
- Invalid input was rejected with a readable message.
- DOM-level tests covered reference loading under a GitHub-style URL subpath,
  input edits, tab switching, mass tables, JSON application, and the asynchronous
  calculation lifecycle using a mocked PyScript transport.
- The full site was not exercised in actual Chrome/Edge/Firefox/Safari here;
  multi-browser visual and PyScript-worker integration verification remains to
  be done. The public hosted link makes that final check convenient.

Read `MODEL_GUIDE.md` for engineering assumptions, equation mapping and units.
The fighter landing-gear equations remain pending verification against their
original source page; other limitations from the desktop teaching model remain.

## Official references

- GitHub publishing source: https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site
- GitHub quick start: https://docs.github.com/en/pages/quickstart
- PyScript configuration: https://docs.pyscript.net/2026.7.3/user-guide/configuration/
- PyScript JavaScript integration: https://docs.pyscript.net/2026.7.3/user-guide/from_javascript/
