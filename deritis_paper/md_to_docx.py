"""Turn manuscript_full.md into the Word file, inserting tables and figures."""

from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

ROOT = Path(Path(__file__).resolve().parent / "deritis_paper")
SRC = ROOT / "manuscript_full.md"
OUT = ROOT / "AST_ALT_sepsis_draft.docx"


def add_para(doc, text, center=False, size=12, bold=False):
    import re
    p = doc.add_paragraph()
    if center:
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pieces = re.split(r"(\*\*[^*]+\*\*)", text) if not bold else [text]
    for piece in pieces:
        if not piece:
            continue
        is_bold = bold or (piece.startswith("**") and piece.endswith("**"))
        run = p.add_run(piece[2:-2] if piece.startswith("**") else piece)
        run.bold = is_bold
        run.font.name = "Times New Roman"
        run.font.size = Pt(size)
    p.paragraph_format.space_after = Pt(8)
    p.paragraph_format.line_spacing = 1.15
    return p


def add_table(doc, headers, rows):
    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.style = "Table Grid"
    for j, h in enumerate(headers):
        table.rows[0].cells[j].text = h
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            table.rows[i + 1].cells[j].text = str(value)
    doc.add_paragraph()


def emit_line(doc, line):
    if line.startswith("Figure ") and ". " in line[:12]:
        number = line.split(".", 1)[0].split()[-1]
        picture = ROOT / "figures" / f"Figure_{number}.png"
        if picture.exists():
            doc.add_picture(str(picture), width=Inches(6.3))
    add_para(doc, line, center=line.startswith("Xiaoye"))


def main():
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(12)
    lines = [raw.strip() for raw in SRC.read_text(encoding="utf-8").splitlines() if raw.strip()]
    sections = []
    front = []
    current = None
    for line in lines:
        if line.startswith("### "):
            current = {"level": 2, "title": line[4:], "body": []}
            sections.append(current)
        elif line.startswith("## "):
            current = {"level": 1, "title": line[3:], "body": []}
            sections.append(current)
        elif line.startswith("# "):
            front.append(("title", line[2:]))
        elif current is None:
            front.append(("text", line))
        else:
            current["body"].append(line)

    body_words = []
    counting = False
    for section in sections:
        if section["title"] == "Introduction":
            counting = True
        if section["title"] == "Acknowledgments":
            counting = False
        if counting:
            body_words.extend(section["body"])
    abstract_words = [line for line in next(s["body"] for s in sections if s["title"] == "Abstract") if not line.lower().startswith("keywords")]
    main_n = len(" ".join(body_words).split())
    abstract_n = len(" ".join(abstract_words).split())

    seen_title = False
    for kind, text in front:
        if kind == "title":
            if not seen_title:
                add_para(doc, text, center=True, size=16, bold=True)
                seen_title = True
            else:
                add_para(doc, text)
            continue
        if text.lower().startswith("short title:"):
            doc.add_heading("Short title", level=1)
            add_para(doc, text.split(":", 1)[1].strip())
            continue
        if text.lower().startswith("article type:"):
            doc.add_heading("Article type", level=1)
            add_para(doc, text.split(":", 1)[1].strip())
            add_para(doc, f"Main-text word count: {main_n:,}; Abstract word count: {abstract_n:,}; Figures: 4; Tables: 6.")
            continue
        add_para(doc, text, center=text.startswith("Xiaoye"))

    figure_section = None
    for section in sections:
        if section["title"] == "Figure legends":
            figure_section = section
            continue
        if section["title"] == "Ethics statement" and section["level"] == 2:
            continue
        title = "Materials and Methods" if section["title"] == "Methods" else section["title"]
        doc.add_heading(title, level=section["level"])
        for line in section["body"]:
            emit_line(doc, line)
        if section["title"] == "Author contributions":
            ethics = next(s for s in sections if s["title"] == "Ethics statement")
            doc.add_heading("Ethics statement", level=1)
            for line in ethics["body"]:
                emit_line(doc, line)

    if figure_section:
        doc.add_heading("Figure legends", level=1)
        for line in figure_section["body"]:
            emit_line(doc, line)

    doc.add_heading("Main table", level=1)
    add_para(doc, "Table 1. MIMIC-IV patients with a usable AST/ALT ratio. Deaths are counted from the enzyme day.")
    add_table(
        doc,
        ["", "Overall", "Q1", "Q2", "Q3", "Q4"],
        [
            ["N", "16616", "4157", "4170", "4136", "4153"],
            ["Age, years", "68.9 (57.4-79.6)", "66.4 (54.7-77.1)", "69.3 (57.9-79.9)", "70.8 (60.0-80.8)", "69.5 (57.9-80.1)"],
            ["Male", "9400 (56.6%)", "2542 (61.2%)", "2371 (56.9%)", "2221 (53.7%)", "2266 (54.6%)"],
            ["SOFA", "6 (4-9)", "5 (3-8)", "6 (4-8)", "6 (4-9)", "7 (4-10)"],
            ["Lactate, mmol/L", "1.9 (1.3-3.1)", "1.6 (1.2-2.6)", "1.8 (1.2-3.0)", "1.9 (1.3-3.3)", "2.2 (1.4-3.7)"],
            ["Bilirubin, mg/dL", "0.6 (0.4-1.2)", "0.6 (0.4-1.3)", "0.6 (0.4-1.1)", "0.6 (0.4-1.1)", "0.7 (0.4-1.3)"],
            ["Creatinine, mg/dL", "1.2 (0.8-2.0)", "1.1 (0.8-1.8)", "1.2 (0.8-1.9)", "1.2 (0.8-2.0)", "1.3 (0.9-2.3)"],
            ["AST, U/L", "40 (23-92)", "35 (21-75)", "33 (21-75)", "38 (23-80)", "60 (33-145)"],
            ["ALT, U/L", "28 (16-65)", "48 (27-111)", "27 (17-62)", "22 (13-47)", "21 (12-46)"],
            ["AST/ALT", "1.4 (1.0-2.0)", "0.8 (0.6-0.9)", "1.2 (1.1-1.3)", "1.7 (1.6-1.8)", "2.7 (2.3-3.4)"],
            ["28-day death", "4249 (25.6%)", "825 (19.8%)", "938 (22.5%)", "1135 (27.4%)", "1351 (32.5%)"],
        ],
    )
    add_para(doc, "Note. n = 16,616 after removal of 5 stays whose enzyme chart time followed the hospital death time. Death is within 28 days of the enzyme measurement. Continuous variables are median (interquartile range).")
    add_para(doc, "Table 2. Analysis cohort and eligible MIMIC-IV stays without a valid AST/ALT pair. Deaths are counted from ICU admission.")
    add_table(
        doc,
        ["", "Usable AST and ALT", "No usable pair"],
        [
            ["N", "16616", "16780"],
            ["Age, years", "68.9 (57.4-79.6)", "70.0 (59.2-80.1)"],
            ["Male", "9400 (56.6%)", "9591 (57.2%)"],
            ["SOFA", "6 (4-9)", "4 (3-6)"],
            ["Lactate, mmol/L", "1.9 (1.3-3.1)", "1.8 (1.3-2.7)"],
            ["Creatinine, mg/dL", "1.2 (0.8-2.0)", "1.0 (0.7-1.5)"],
            ["Cardiac surgical ICU", "1294 (7.8%)", "4760 (28.4%)"],
            ["28-day death from ICU admission", "4239 (25.5%)", "2285 (13.6%)"],
        ],
    )
    add_para(doc, "Note. The first column is the analysis cohort of 16,616, after removal of 5 stays whose enzyme chart time followed hospital death. The second column is eligible stays with no valid enzyme pair. Death in both columns is counted from ICU admission.")
    add_para(doc, "Table 3. Adjusted associations.")
    add_table(
        doc,
        ["Contrast", "Estimate", "95% CI", "P"],
        [
                    ["Q2 vs Q1, imputed", "HR 1.04", "0.95-1.14", "0.41"],
                    ["Q3 vs Q1, imputed", "HR 1.17", "1.07-1.29", "<0.001"],
                    ["Q4 vs Q1, imputed", "HR 1.43", "1.31-1.57", "<0.001"],
            ["Q4 vs Q1, complete case", "HR 1.36", "1.23-1.51", "<0.001"],
                    ["Q4 vs Q1, days 0-7", "HR 1.81", "1.60-2.05", "<0.001"],
                    ["Q4 vs Q1, days 7-28", "HR 1.13", "0.99-1.28", "0.065"],
                    ["Days 7-28 relative to days 0-7", "HR 0.62", "0.52-0.74", "<0.001"],
            ["Log AST per SD, imputed", "HR 1.42", "1.33-1.52", "<0.001"],
            ["Log ALT per SD, imputed", "HR 0.77", "0.72-0.82", "<0.001"],
            ["Q4 vs Q1, imputed, expanded covariates", "HR 1.33", "1.22-1.46", "<0.001"],
            ["Q4 vs Q1, excluding alcohol, hepatitis, fatty liver", "HR 1.40", "1.25-1.56", "<0.001"],
            ["eICU Q2 vs Q1, imputed, clustered SE", "OR 1.15", "0.99-1.35", "0.068"],
            ["eICU Q3 vs Q1, imputed, clustered SE", "OR 1.30", "1.14-1.49", "<0.001"],
            ["eICU Q4 vs Q1, imputed, clustered SE", "OR 1.95", "1.69-2.26", "<0.001"],
        ],
    )
    add_para(doc, "Note. Quartile rows for 28-day death, and the log AST and log ALT rows, are multiply imputed on 16,616 patients. The day 0-7, day 7-28, and relative rows are from one time-varying model on the same 30 imputed datasets. The complete-case quartile row and the complete-case joint enzyme estimates use 10,890 patients; those joint estimates are also in Table 5. The eICU rows use all 12,374 patients and hospital-clustered standard errors. The complete-case eICU estimate is in Table 5.")
    add_para(doc, "Table 4. Interaction hazard ratios in one sample: 5,423 complete cases from the lowest and highest quartiles, with 1,629 deaths. Every row uses this sample.")
    add_table(
        doc,
        ["Contrast", "HR", "95% CI"],
        [
            ["Men versus women", "1.03", "0.84-1.27"],
            ["Age >= 65 versus younger", "0.80", "0.64-1.00"],
            ["Lactate >= 4 versus lower", "1.09", "0.85-1.41"],
            ["Cardiac surgical ICU versus other units", "0.61", "0.36-1.03"],
        ],
    )
    add_para(doc, "Note. N = 5,423 and deaths = 1,629 are the size of the Q1-plus-Q4 complete-case analysis. A value below 1 means a smaller association in the named group.")
    add_para(doc, "Table 5. Sensitivity analyses for timing, selection, enzyme bounds, and hospital death.")
    add_table(
        doc,
        ["Analysis", "Contrast", "Estimate", "95% CI"],
        [
            ["Alive at 24 hours, days 0-7", "Q4 vs Q1", "HR 1.82", "1.59-2.08"],
            ["Measurement-day deaths removed, days 0-7", "Q4 vs Q1", "HR 1.81", "1.59-2.06"],
            ["Measurement-day deaths removed, days 7-28", "Q4 vs Q1", "HR 1.15", "1.01-1.30"],
            ["Measurement-day deaths removed, later relative to early", "Q4", "HR 0.63", "0.53-0.76"],
            ["Alive at 24 hours, days 7-28", "Q4 vs Q1", "HR 1.17", "1.03-1.32"],
            ["Days 0-3, imputed", "Q4 vs Q1", "HR 2.04", "1.72-2.41"],
            ["After day 3, imputed", "Q4 vs Q1", "HR 1.25", "1.12-1.38"],
            ["Days 0-14, imputed", "Q4 vs Q1", "HR 1.62", "1.46-1.79"],
            ["After day 14, imputed", "Q4 vs Q1", "HR 0.99", "0.83-1.19"],
            ["Inverse-probability weighted", "Q4 vs Q1", "HR 1.46", "1.32-1.61"],
            ["AST/ALT bounds 1-10,000, complete case", "Q4 vs Q1", "HR 1.33", "1.20-1.47"],
            ["AST/ALT bounds 10-1,000, complete case", "Q4 vs Q1", "HR 1.41", "1.27-1.58"],
            ["MIMIC-IV hospital death, imputed", "Q4 vs Q1", "OR 1.49", "1.32-1.67"],
            ["Log AST per SD, complete case", "AST per SD", "HR 1.38", "1.28-1.48"],
            ["Log ALT per SD, complete case", "ALT per SD", "HR 0.81", "0.75-0.87"],
            ["Log ratio per SD, imputed", "Log AST minus log ALT", "HR 1.16", "1.13-1.20"],
            ["Log mean level per SD, imputed", "Mean of log AST and log ALT", "HR 1.09", "1.06-1.13"],
            ["eICU complete case, hospital-clustered SE", "Q4 vs Q1", "OR 1.96", "1.63-2.36"],
        ],
    )
    add_para(doc, "Note. Quartile indicators were predictors in the imputation model. The rows labeled “Alive at 24 hours” exclude patients who died within 24 hours after ICU admission. Follow-up is still counted from the day of enzyme measurement. The measurement-day rows remove the 269 deaths on the calendar day of the enzyme pair from the 16,616-patient cohort and refit the day-7 split. The ratio row is the later hazard divided by the first-week hazard. The log ratio and log mean rows are the joint enzyme model rewritten in that form; the mean-level term is the test of whether the ratio alone is sufficient. The day-3 and day-14 rows use that same time-varying model. Inverse-probability weights used age, sex, SOFA, alcohol use disorder, and care unit, and were applied to all 16,616 analyzed stays; that model does not include lactate. Enzyme-bound rows keep the original quartile cutpoints. The hospital-death odds ratio uses the same covariates as the 28-day hazard model. The log AST and log ALT rows are the complete-case joint model (10,890 patients). The eICU row is the complete-case model with standard errors clustered by hospital (157 hospitals).")
    add_para(doc, "Table 6. Missing covariates by AST/ALT quartile in MIMIC-IV.")
    add_table(
        doc,
        ["Variable", "Overall", "Q1", "Q2", "Q3", "Q4"],
        [
            ["Lactate missing", "31.9%", "40.0%", "33.1%", "29.9%", "24.8%"],
            ["Bilirubin missing", "3.7%", "3.9%", "4.1%", "3.6%", "3.3%"],
            ["Creatinine missing", "0.1%", "0.1%", "0.1%", "0.1%", "0.1%"],
            ["SOFA missing", "0", "0", "0", "0", "0"],
            ["Albumin missing", "42.3%", "42.8%", "42.1%", "42.3%", "42.0%"],
        ],
    )
    add_para(doc, "Note. Percent of the 16,616 analyzed stays. Lactate is the covariate that removes the most patients from a complete-case model, and it is missing more often in the lower quartiles. Alcohol use disorder, age, and sex had no missing values. Multiple imputation assumes that, given the variables in the imputation model, the missing values do not depend on the unobserved number.")
    doc.save(OUT)
    text = "\n".join(p.text for p in doc.paragraphs)
    print("paras", len(doc.paragraphs), "tables", len(doc.tables), "pics", len(doc.inline_shapes))
    print("refs", text.count("doi:"))


if __name__ == "__main__":
    main()
