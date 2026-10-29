# Public documentation

The site at https://diptavo.github.io/PASSAGE/ is built from the HTML fragments
in `pages/`, a shared template and local CSS/JavaScript assets. It has no runtime
dependencies or external fonts. Statistical details remain in `docs/`; the
website provides an accessible methods overview and a practical software guide.

Build from the repository root with Python 3.9 or later:

```bash
python3 site/build.py
python3 -m http.server 8000 --directory site/_build --bind 127.0.0.1
```

Open http://127.0.0.1:8000/ to preview. The build checks internal links, anchors,
duplicate IDs and page headings. Relative links support the `/PASSAGE/` project
prefix on GitHub Pages. Generated output is ignored by Git and excluded from
the R package.

The Documentation workflow validates pull requests and publishes changes on
`main`. Repository Settings → Pages must use GitHub Actions as the source.
Only `site/_build` is uploaded; local analysis files and other documentation
directories are not part of the deployment.
