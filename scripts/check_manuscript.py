#!/usr/bin/env python
"""Check the drafted manuscript against the style rules agreed for this paper.

Six rules came out of drafting Sections 2.1 and 2.2 against the author's own papers. Four of
them are mechanical and are enforced here. The other two are matters of judgement and are
only reported on, because a script that pretends to check them would be worse than one that
says it cannot.

  1. A drafted section carries no \\TODO and no \\NOTE.                        [reported]
  2. Paragraphs are not left short.                                           [reported]
  3. No colons in prose or in captions.                                       [enforced]
  4. No "rung" or "ladder" before the subsection that defines them.           [enforced]
  5. Captions identify what a panel shows; numbers belong in the prose.       [reported]
  6. Every figure panel is referenced explicitly in the prose.                [enforced]
  7. Procedure in Methods is written in the past tense.                       [reported]
  8. Sentences stay between 12 and 45 words.                                  [reported]

Rule 6 exists because all three figures in Sections 2.1 and 2.2 were first drafted with the
figure referenced only as a whole. The content of every panel was discussed, so nothing was
missing from the argument, and a reader still had no way to tell which panel to look at. That
is the kind of defect that survives proofreading, which is why it is checked by a script.

Exit status is non-zero if any enforced rule fails, so this can run before a commit.

Usage:
    python scripts/check_manuscript.py
    python scripts/check_manuscript.py --sections 2.1 2.2
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SEC = ROOT/"manuscript"/"sections"
RUNG_DEFINED_IN = "increasing computational"   # the subsection that introduces them
# Rule 4 covers the vocabulary of the descriptor ladder, which is this project's own.
# "rung" and "ladder" are the obvious members; "molecular-dynamics descriptor" joins
# them because it presupposes that the sets have been introduced and ordered, which
# only that subsection does. The simulations themselves may be named earlier, since
# the opening roadmap has to say what the section contains.
EARLY_TERMS = r"\b(rungs?|ladder)\b|molecular[- ]dynamics descriptor"
SHORT_PARA_WORDS = 60
SENT_MIN, SENT_MAX = 12, 45

# Splitting on ". " alone breaks at "i.e." and at siunitx remnants, which invents
# short sentences that are not there. These are masked before the split.
ABBREV = (r"i\.e\.", r"e\.g\.", r"cf\.", r"et al\.", r"vs\.", r"Fig\.", r"Eq\.",
          r"approx\.", r"\b[A-Z]\.")


def sentences(text: str) -> list:
    t = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})*", " ", text)
    t = re.sub(r"[${}~^_\\]", " ", t)
    t = re.sub(r"(\d)\.(\d)", r"\1<DOT>\2", t)
    for a in ABBREV:
        t = re.sub(a, lambda m: m.group(0).replace(".", "<DOT>"), t)
    t = re.sub(r"\s+", " ", t)
    out = [s.replace("<DOT>", ".").strip() for s in re.split(r"(?<=[.!?]) +", t)]
    return [s for s in out if len(s.split()) > 1]

# Rule 7. Methods describes what was done, so a present-tense passive there is
# usually a slip. Three kinds of present tense are correct and are allowed. A
# pointer to the paper's own apparatus ("is shown in Figure~2"), a timeless fact
# about the protein or the chemistry ("Trp86 is located at the catalytic anionic
# site"), and a definition that still holds ("the descriptors are defined in terms
# of residue contacts"). The allowance is by verb rather than by context, so the
# check reports rather than enforces and a few false positives are expected.
PRESENT_PASSIVE = re.compile(r"\b(?:is|are)\s+(?:\w+ly\s+)?(\w+ed|shown|drawn|built|"
                             r"set|kept|taken|made|written|held|run|known)\b")
TENSE_ALLOWED = {
    "shown", "described", "compared", "drawn", "listed", "summarised", "summarized",
    "located", "situated", "designed", "reported", "detected", "marked", "known",
    "defined", "resolved", "referred", "based", "expressed", "formed", "composed",
    "related", "concerned", "intended", "regarded", "considered", "called", "termed",
}


def figure_blocks(text: str, env: str = "figure"):
    """(label, caption, full block) for every float environment of this kind."""
    for m in re.finditer(r"\\begin\{" + env + r"\}(.*?)\\end\{" + env + r"\}", text, re.S):
        b = m.group(1)
        lab = re.search(r"\\label\{((?:fig|tab):[a-zA-Z]+)\}", b)
        i = b.find(r"\caption{")
        if i < 0 or not lab:
            continue
        j = i + len(r"\caption{")
        depth, k = 1, j
        while k < len(b) and depth:
            depth += (b[k] == "{") - (b[k] == "}")
            k += 1
        yield lab.group(1), " ".join(b[j:k-1].split()), m.group(0)


def drop_macro(text: str, name: str) -> str:
    """Remove every \\name{...}, counting braces.

    A regex with one level of nesting baked in was the first attempt and it silently
    failed on \\TODO{... $\\sigma^2\\_{\\mathrm{lig}}$ ...}, which is two levels deep, so
    that TODO was read as prose by every rule that works on prose.
    """
    out, i, tag = [], 0, "\\" + name + "{"
    while True:
        j = text.find(tag, i)
        if j < 0:
            out.append(text[i:])
            return "".join(out)
        out.append(text[i:j])
        k, depth = j + len(tag), 1
        while k < len(text) and depth:
            depth += (text[k] == "{") - (text[k] == "}")
            k += 1
        out.append(" ")
        i = k


def prose_of(text: str) -> str:
    """Everything that is not a float environment, a comment or a drafting macro."""
    t = re.sub(r"\\begin\{(figure|table)\}.*?\\end\{\1\}", " ", text, flags=re.S)
    t = re.sub(r"(?m)^\s*%.*$", " ", t)
    for name in ("TODO", "NOTE"):
        t = drop_macro(t, name)
    return t


def strip_refs(t: str) -> str:
    return re.sub(r"\\(ref|label|cite[tp]?|citep|citet)\{[^}]*\}", "", t)


def strip_headings(t: str) -> str:
    """Section and paragraph headings are titles, not sentences, and the colon rule is
    about sentences. Removing them here keeps the check from reporting a style choice
    that was never in question."""
    return re.sub(r"\\(sub)*(section|paragraph)\*?\{[^}]*\}", " ", t)


def subsection_prose(text: str) -> dict:
    """{subsection title: its prose}, with float environments and macros removed."""
    out, parts = {}, re.split(r"(\\subsection\*?\{[^}]*\})", text)
    cur = "(preamble)"
    for chunk in parts:
        m = re.match(r"\\subsection\*?\{([^}]*)\}", chunk)
        if m:
            cur = m.group(1)
            out.setdefault(cur, "")
            continue
        out[cur] = out.get(cur, "") + " " + prose_of(chunk)
    return out


def drafted_subsections(text: str) -> dict:
    """{subsection title: has real prose} -- a figure in a subsection that is still
    scaffolding cannot yet reference its panels, and flagging it is noise."""
    return {k: len(v.split()) > 60 for k, v in subsection_prose(text).items()}


def panels_cited(label: str, prose: str, sub_prose: str) -> set:
    """Panels a reader is actually pointed at.

    Two forms count. `Figure~\\ref{fig:x}(a)` anywhere, and a bare `panel~(a)` inside
    the subsection that already names the figure. The second form was not accepted at
    first, which failed two honest paragraphs: naming the figure once and then writing
    "the temperature in panel~(a)" reads better than repeating "Figure~5(a)" four times
    in four sentences, and it points the reader just as precisely.
    """
    cited = set(re.findall(r"\\ref\{" + re.escape(label) + r"\}\(([a-dA-D])\)", prose))
    if re.search(r"\\ref\{" + re.escape(label) + r"\}", sub_prose):
        # case-insensitive, since a sentence may open with "Panel~(d) shows ..."
        cited |= set(re.findall(r"panels?~?\s*\(([a-dA-D])\)", sub_prose, re.I))
    return cited


def subsection_of(text: str, pos: int) -> str:
    ms = list(re.finditer(r"\\subsection\*?\{([^}]*)\}", text[:pos]))
    return ms[-1].group(1) if ms else "(preamble)"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--files", nargs="*", default=None)
    a = ap.parse_args()
    files = [Path(f) for f in a.files] if a.files else sorted(SEC.glob("*.tex"))

    failures, notes = [], []
    for f in files:
        text = f.read_text()
        prose = prose_of(text)
        drafted = len(prose.split()) > 120 and r"\TODO" not in text

        # rule 1 and 2 -- reported, not enforced
        n_todo = text.count(r"\TODO")
        n_note = text.count(r"\NOTE")
        if drafted and (n_todo or n_note):
            notes.append(f"{f.name}: drafted but still carries {n_todo} TODO, {n_note} NOTE")
        for para in [p for p in re.split(r"\n\s*\n", prose) if p.strip()]:
            w = len(strip_refs(para).split())
            if 0 < w < SHORT_PARA_WORDS and not para.strip().startswith("\\"):
                notes.append(f"{f.name}: short paragraph ({w} words) -- "
                             f"{' '.join(para.split())[:60]}...")

        # rule 3 -- colons anywhere a reader reads
        for label, cap, _blk in figure_blocks(text):
            if ":" in strip_refs(cap):
                failures.append(f"{f.name}: colon in the caption of {label}")
        for line in strip_headings(strip_refs(prose)).split("\n"):
            if ":" in line:
                failures.append(f"{f.name}: colon in prose -- {line.strip()[:70]}")

        # rule 4 -- vocabulary introduced before it is defined
        cut = text.find(RUNG_DEFINED_IN)
        if cut > 0:
            head = strip_refs(prose_of(text[:cut]))
            for m in re.finditer(EARLY_TERMS, head, re.I):
                ctx = " ".join(head[max(0, m.start()-45):m.start()+25].split())
                failures.append(f"{f.name}: '{m.group(0)}' used before "
                                f"'{RUNG_DEFINED_IN}' -- ...{ctx}...")

        # rule 5 -- caption length, reported. A table caption is held to a much
        # shorter limit than a figure's. A table caption should name what the
        # columns hold, so that the table stands on its own, and stop there; the
        # qualifications and the interpretation belong in the prose. Four words is
        # too few and a hundred is prose in the wrong place.
        for env, limit in (("figure", 120), ("table", 35)):
            for label, cap, _blk in figure_blocks(text, env):
                n = len(cap.split())
                if n > limit:
                    notes.append(f"{f.name}: caption of {label} is {n} words, over the "
                                 f"{limit} allowed for a {env}; move it into the prose")

        # rule 7 -- Methods in the past tense, reported
        if "methods" in f.name:
            for line in strip_headings(strip_refs(prose)).split("\n"):
                for m in PRESENT_PASSIVE.finditer(line):
                    if m.group(1) in TENSE_ALLOWED:
                        continue
                    notes.append(f"{f.name}: present-tense passive in Methods -- "
                                 f"...{line.strip()[max(0, m.start()-35):m.end()+15]}...")

        # rule 8 -- sentence length, reported. The reference papers run 22 to 29
        # words per sentence on average; the risk here is variance rather than the
        # mean, i.e. a four-word sentence next to an eighty-word one.
        for s in sentences(strip_headings(strip_refs(prose))):
            n = len(s.split())
            if n < SENT_MIN:
                notes.append(f"{f.name}: short sentence ({n} words) -- {s}")
            elif n > SENT_MAX:
                notes.append(f"{f.name}: long sentence ({n} words) -- {s[:70]}...")

        # rule 6 -- every panel referenced, but only once the subsection has prose
        subs = subsection_prose(text)
        done = {k: len(v.split()) > 60 for k, v in subs.items()}
        for label, cap, blk in figure_blocks(text):
            panels = sorted(set(re.findall(r"\((?:~)?([a-dA-D])\)", cap)))
            if not panels:
                continue
            sub = subsection_of(text, text.find(blk))
            cited = panels_cited(label, prose, subs.get(sub, ""))
            missing = [p for p in panels if p not in cited]
            if not missing:
                continue
            msg = f"{f.name}: {label} panel(s) {','.join(missing)} never referenced in prose"
            (failures if done.get(sub) else notes).append(
                msg + ("" if done.get(sub) else f"  [{sub} not drafted yet]"))

    for n in notes:
        print(f"  note     {n}")
    for x in failures:
        print(f"  FAIL     {x}")
    print(f"\n  {len(failures)} failure(s), {len(notes)} note(s)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
