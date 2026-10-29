#!/usr/bin/env python3
"""Build the public documentation using only the Python standard library."""
from html import escape
from html.parser import HTMLParser
from pathlib import Path
import re
import shutil
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "_build"
PAGES = {
    "index": ("PASSAGE", "Spatial pathway analysis, from the scientific question to an interpretable result."),
    "methods": ("Methods", "The questions, estimation choices and inference behind PASSAGE."),
    "software": ("Software guide", "Install PASSAGE, prepare an input bundle, estimate pathway activity and run both testing workflows."),
    "example": ("Example analysis", "Analyze three Hallmark pathways in the prepared KC2 kidney cancer sample."),
}


class Document(HTMLParser):
    def __init__(self, source):
        super().__init__(convert_charrefs=True)
        self.ids, self.links, self.headings = set(), [], []
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            if attrs["id"] in self.ids:
                raise ValueError("Duplicate id: " + attrs["id"])
            self.ids.add(attrs["id"])
        for key in ("href", "src"):
            if key in attrs:
                self.links.append(attrs[key])
        if tag == "h1":
            self.headings.append(tag)


def build():
    OUTPUT.mkdir(exist_ok=True)
    shutil.copytree(ROOT / "assets", OUTPUT / "assets", dirs_exist_ok=True)
    template = (ROOT / "template.html").read_text()
    for slug, (title, description) in PAGES.items():
        prefix = "./" if slug == "index" else "../"
        body = (ROOT / "pages" / (slug + ".html")).read_text()
        headings = re.findall(r'<h2 id="([^"]+)">(.+?)</h2>', body)
        toc = "".join('<a href="#{}">{}</a>'.format(anchor, text) for anchor, text in headings)
        nav = "".join('<a {}href="{}{}">{}</a>'.format(
            'aria-current="page" ' if slug == key else "", prefix, key + "/", label
        ) for key, label in (("methods", "Methods"), ("software", "Software guide"), ("example", "Example analysis")))
        values = {"title": title, "description": escape(description), "body": body,
                  "nav": nav, "toc": toc, "prefix": prefix,
                  "route": "" if slug == "index" else slug + "/"}
        rendered = template
        for key, value in values.items():
            rendered = rendered.replace("{{" + key + "}}", value)
        # Body fragments also use the shared root-relative-to-page prefix.
        rendered = rendered.replace("{{prefix}}", prefix)
        if "{{" in rendered:
            raise ValueError("Unresolved template field in " + slug)
        destination = OUTPUT if slug == "index" else OUTPUT / slug
        destination.mkdir(exist_ok=True)
        (destination / "index.html").write_text(rendered)
    (OUTPUT / ".nojekyll").write_text("")
    validate()


def validate():
    documents = {p.resolve(): Document(p.read_text()) for p in OUTPUT.rglob("*.html")}
    checked = 0
    for path, doc in documents.items():
        if len(doc.headings) != 1:
            raise ValueError(str(path) + " must contain exactly one h1")
        for link in doc.links:
            url = urlsplit(link)
            if url.scheme or url.netloc:
                continue
            if url.path.startswith("/"):
                raise ValueError("Use project-safe relative links: " + link)
            target = (path.parent / unquote(url.path)).resolve() if url.path else path
            if target.is_dir():
                target = target / "index.html"
            if not target.is_relative_to(OUTPUT.resolve()) or not target.exists():
                raise ValueError("Broken link in {}: {}".format(path.name, link))
            if url.fragment and unquote(url.fragment) not in documents.get(target, doc).ids:
                raise ValueError("Missing anchor: " + link)
            checked += 1
    print("Built {} pages; verified {} internal links and assets.".format(len(documents), checked))


if __name__ == "__main__":
    build()
