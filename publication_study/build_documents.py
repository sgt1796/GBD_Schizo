from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


BLUE = "1F4E78"
LIGHT_BLUE = "D9EAF7"
LIGHT_ORANGE = "FCE4D6"
GRAY = "666666"
GREEN = "E2F0D9"


def set_cell_shading(cell, fill: str):
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_margins(cell, top=80, start=90, bottom=80, end=90):
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for margin, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{margin}"))
        if node is None:
            node = OxmlElement(f"w:{margin}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value)); node.set(qn("w:type"), "dxa")


def repeat_table_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def set_repeat_table_header(row):
    repeat_table_header(row)


def prevent_row_split(row):
    tr_pr = row._tr.get_or_add_trPr()
    cant_split = OxmlElement("w:cantSplit")
    tr_pr.append(cant_split)


def add_page_number(paragraph):
    paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    run = paragraph.add_run("Page ")
    fld = OxmlElement("w:fldSimple"); fld.set(qn("w:instr"), "PAGE")
    run._r.addnext(fld)


def set_image_alt(inline_shape, text: str):
    doc_pr = inline_shape._inline.docPr
    doc_pr.set("descr", text)
    doc_pr.set("title", text[:80])


def setup_document(title: str, subtitle: str | None = None, provisional: bool = False, show_blinded: bool = True) -> Document:
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Inches(0.7); section.bottom_margin = Inches(0.7)
    section.left_margin = Inches(0.8); section.right_margin = Inches(0.8)
    styles = doc.styles
    normal = styles["Normal"]
    normal.font.name = "Arial"; normal.font.size = Pt(10.5)
    normal.paragraph_format.space_after = Pt(6); normal.paragraph_format.line_spacing = 1.12
    for name, size, color in (("Title", 20, BLUE), ("Heading 1", 14, BLUE), ("Heading 2", 11.5, BLUE), ("Heading 3", 10.5, GRAY)):
        style = styles[name]; style.font.name = "Arial"; style.font.size = Pt(size); style.font.color.rgb = RGBColor.from_string(color)
        style.font.bold = True; style.paragraph_format.space_before = Pt(10); style.paragraph_format.space_after = Pt(5)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(title); r.bold = True; r.font.name = "Arial"; r.font.size = Pt(20); r.font.color.rgb = RGBColor.from_string(BLUE)
    if subtitle:
        p2 = doc.add_paragraph(); p2.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r2 = p2.add_run(subtitle); r2.italic = True; r2.font.size = Pt(11); r2.font.color.rgb = RGBColor.from_string(GRAY)
    if show_blinded:
        p3 = doc.add_paragraph(); p3.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p3.add_run("Blinded manuscript").bold = True
    if provisional:
        q = doc.add_paragraph(); q.alignment = WD_ALIGN_PARAGRAPH.CENTER
        q.paragraph_format.left_indent=Inches(.05); q.paragraph_format.right_indent=Inches(.05)
        q.paragraph_format.space_before=Pt(4); q.paragraph_format.space_after=Pt(10)
        q_pr=q._p.get_or_add_pPr(); shd=OxmlElement("w:shd"); shd.set(qn("w:fill"),LIGHT_ORANGE); q_pr.append(shd)
        run = q.add_run("PROVISIONAL ANALYTICAL DRAFT - NOT FOR SUBMISSION\nOfficial GBD 2023 population and official NCI Joinpoint outputs remain external submission gates.")
        run.bold = True; run.font.color.rgb = RGBColor(156, 87, 0)
    footer = section.footer.paragraphs[0]
    add_page_number(footer)
    return doc


def add_heading(doc: Document, text: str, level: int = 1):
    return doc.add_heading(text, level=level)


def add_paragraph(doc: Document, text: str, bold_lead: str | None = None):
    p = doc.add_paragraph()
    if bold_lead and text.startswith(bold_lead):
        p.add_run(bold_lead).bold = True
        p.add_run(text[len(bold_lead):])
    else:
        p.add_run(text)
    return p


def add_bullets(doc: Document, items: list[str]):
    for item in items:
        doc.add_paragraph(item, style="List Bullet")


def add_dataframe_table(doc: Document, frame: pd.DataFrame, columns: list[tuple[str, str]], font_size=8.2, widths=None):
    table = doc.add_table(rows=1, cols=len(columns)); table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = "Table Grid"
    hdr = table.rows[0]; repeat_table_header(hdr)
    for i, (_, label) in enumerate(columns):
        cell = hdr.cells[i]; cell.text = label; set_cell_shading(cell, BLUE); set_cell_margins(cell)
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        for run in cell.paragraphs[0].runs:
            run.bold = True; run.font.color.rgb = RGBColor(255,255,255); run.font.size = Pt(font_size)
        cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
    for _, row in frame.iterrows():
        cells = table.add_row().cells
        prevent_row_split(table.rows[-1])
        for i, (key, _) in enumerate(columns):
            value = row.get(key, "")
            cells[i].text = "" if pd.isna(value) else str(value)
            cells[i].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER; set_cell_margins(cells[i])
            for run in cells[i].paragraphs[0].runs: run.font.size = Pt(font_size)
            cells[i].paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.LEFT if i == 0 else WD_ALIGN_PARAGRAPH.CENTER
            if widths: cells[i].width = Inches(widths[i])
    doc.add_paragraph()
    return table


def add_figure(doc: Document, path: Path, caption: str, alt: str, width=6.7):
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    shape = p.add_run().add_picture(str(path), width=Inches(width)); set_image_alt(shape, alt)
    cap = doc.add_paragraph(); cap.alignment = WD_ALIGN_PARAGRAPH.LEFT
    lead, rest = caption.split(". ", 1)
    cap.add_run(lead + ". ").bold = True; cap.add_run(rest)
    cap.paragraph_format.space_after = Pt(10)


def fmt(value, digits=1):
    if pd.isna(value): return "NA"
    return f"{float(value):,.{digits}f}"


def endpoint_lookup(endpoints: pd.DataFrame, loc: str, sex: str, outcome: str, metric: str):
    return endpoints[(endpoints.location_name==loc)&(endpoints.sex_name==sex)&(endpoints.measure_name==outcome)&(endpoints.metric_name==metric)].iloc[0]


def compact_endpoint_table(endpoints: pd.DataFrame, segmented: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for loc in ("China","United States of America"):
        for sex in ("Female","Male"):
            for outcome in ("Incidence","Prevalence","DALYs"):
                count=endpoint_lookup(endpoints,loc,sex,outcome,"All-age count")
                rate=endpoint_lookup(endpoints,loc,sex,outcome,"Age-standardized rate per 100,000")
                seg=segmented[(segmented.location_name==loc)&(segmented.sex_name==sex)&(segmented.measure_name==outcome)].iloc[0]
                rows.append({"Location":loc.replace("United States of America","United States"),"Sex":sex,"Outcome":outcome,
                             "Count 2023 (95% UI)":f"{fmt(count.value_2023,0)} ({fmt(count.lower_2023,0)}-{fmt(count.upper_2023,0)})",
                             "Count change, %":fmt(count.percent_change_point_estimate,1),
                             "ASR 1990":fmt(rate.value_1990,2),"ASR 2023 (95% UI)":f"{fmt(rate.value_2023,2)} ({fmt(rate.lower_2023,2)}-{fmt(rate.upper_2023,2)})",
                             "AAPC, %/year":f"{fmt(seg.aapc,3)} ({fmt(seg.aapc_lower_model_ci,3)} to {fmt(seg.aapc_upper_model_ci,3)})"})
    return pd.DataFrame(rows)


def build_manuscript(analysis: Path, out_dir: Path, meta: dict, tables: dict[str,pd.DataFrame]):
    provisional=not meta["submission_ready"]
    doc=setup_document("Divergent sex-specific trends and demographic drivers of schizophrenia burden in China and the United States, 1990-2023",
                       "A comparative analysis of Global Burden of Disease 2023 estimates",provisional)
    endpoints=tables["endpoint_summary"]; segmented=tables["segmented_summary"]; pairs=tables["pairwise_parallelism"]; decomp=tables["decomposition"]

    add_heading(doc,"Abstract")
    china_f=endpoint_lookup(endpoints,"China","Female","DALYs","Age-standardized rate per 100,000")
    china_m=endpoint_lookup(endpoints,"China","Male","DALYs","Age-standardized rate per 100,000")
    us_f=endpoint_lookup(endpoints,"United States of America","Female","DALYs","Age-standardized rate per 100,000")
    us_m=endpoint_lookup(endpoints,"United States of America","Male","DALYs","Age-standardized rate per 100,000")
    significant=int((pairs.q_value_bh_within_family<0.05).sum())
    add_paragraph(doc,"Background: Schizophrenia produces substantial lifelong disability, but comparisons of national trends can be obscured by population growth, population ageing, and changes in age-specific rates. We compared sex-specific schizophrenia burden in China and the United States, two populous settings with contrasting demographic trajectories.")
    add_paragraph(doc,"Methods: We analyzed GBD 2023 incidence, prevalence, and disability-adjusted life-year (DALY) estimates for 1990-2023. We summarized all-age counts and age-standardized rates (ASRs), characterized nonlinear trends with permutation-selected segmented log-linear regression, controlled pairwise comparison multiplicity using the Benjamini-Hochberg procedure, and decomposed adult count changes into population growth, population ageing, and age-specific rate change using permutation-averaged Shapley replacement. Incidence age-period-cohort estimable functions were secondary.")
    add_paragraph(doc,f"Results: From 1990 to 2023, all-age counts increased for every outcome, country, and sex stratum. DALY ASRs in China changed from {fmt(china_f.value_1990,1)} to {fmt(china_f.value_2023,1)} per 100,000 among females and {fmt(china_m.value_1990,1)} to {fmt(china_m.value_2023,1)} among males. Corresponding U.S. rates declined from {fmt(us_f.value_1990,1)} to {fmt(us_f.value_2023,1)} among females and {fmt(us_m.value_1990,1)} to {fmt(us_m.value_2023,1)} among males. After false-discovery-rate control, {significant} of 12 prespecified trajectory comparisons rejected parallelism. Population growth was the principal positive contributor to rising counts in most strata, while age-specific rate change separated the near-stable Chinese rates from declining U.S. rates.")
    add_paragraph(doc,"Conclusions: Rising schizophrenia counts in both countries did not imply worsening age-standardized burden. The China-U.S. divergence was primarily a rate-trajectory and demographic-composition phenomenon, supporting sex- and age-responsive service planning while avoiding causal attribution to health-system differences.")
    add_paragraph(doc,"Keywords: schizophrenia; Global Burden of Disease; China; United States; segmented regression; decomposition; demographic change; sex differences")

    add_heading(doc,"Introduction")
    add_paragraph(doc,"Schizophrenia is a severe mental disorder associated with persistent functional impairment, premature mortality, family burden, and high demand for health and social services. In the Global Burden of Disease (GBD) framework, its direct disease burden is expressed almost entirely through years lived with disability because schizophrenia is not assigned a direct fatal cause component. Consequently, prevalence and disability trends are central to interpreting service need.")
    add_paragraph(doc,"Global analyses consistently show that numbers of people living with schizophrenia and associated DALYs have risen while age-standardized rates have changed much less. This apparent paradox reflects population growth and changing age structures rather than a single epidemiologic process. China and the United States provide an informative comparison: each contains a large population, but their demographic histories and estimated schizophrenia rate trajectories differ.")
    add_paragraph(doc,"Recent GBD publications have already described global schizophrenia patterns, China-specific trends, joinpoint models, decomposition, and forecasts. A further generic trend-and-projection analysis would add little. What remains useful is a prespecified, sex-stratified comparison that formally assesses trajectory non-parallelism and links count changes to demographic and rate components without treating ecological contrasts as causal.")
    add_paragraph(doc,"We therefore examined incidence, prevalence, and DALY burden in China and the United States from 1990 through 2023. Our objectives were to quantify endpoint and nonlinear trends, test country and sex trajectory differences, decompose adult count changes, and determine whether secondary age-period-cohort estimable summaries supported the primary findings.")

    add_heading(doc,"Methods")
    add_heading(doc,"Study design and data source",2)
    add_paragraph(doc,"We conducted a comparative secondary analysis of GBD 2023 modeled health estimates. The analytic population comprised females and males in China and the United States from 1990 through 2023. The burden extract contained annual point estimates and 95% uncertainty intervals (UIs) for incidence, prevalence, years lived with disability (YLDs), and DALYs by age, sex, country, year, and metric. The source files contained aggregated, non-identifiable population estimates; institutional review board review and informed consent were therefore not applicable.")
    add_heading(doc,"Outcomes and prespecified exclusions",2)
    add_paragraph(doc,"Primary outcomes were incidence, prevalence, and DALYs. Trend analyses used all-age numbers and age-standardized rates per 100,000. Decomposition used age-specific rates and populations from ages 15-19 through 70+ years. We excluded the GBD Percent metric because its denominator and age basis differed by outcome. Probability-of-death and available risk-factor extracts were excluded because they did not represent schizophrenia-specific causal attribution. YLD and DALY estimates were audited for numerical identity and YLDs were not duplicated in results.")
    add_heading(doc,"Data quality and population denominators",2)
    source_text="official GBD 2023 population estimates" if meta["population_status"]=="official_GBD_2023" else "a provisional population proxy reconstructed from matched GBD count-rate pairs"
    add_paragraph(doc,f"We audited dimensional uniqueness, completeness across 34 years, missingness, positivity, uncertainty-bound ordering, and age-bin consistency. Decomposition currently uses {source_text}. A production build requires population estimates from the same GBD 2023 release. We checked age-specific reconstruction by multiplying population by rate and dividing by 100,000.")
    add_heading(doc,"Trend and comparison analyses",2)
    trend_method="official NCI Joinpoint 6.0.1" if meta["submission_ready"] else "an independently implemented permutation-selected segmented regression that must not be described as NCI Joinpoint"
    add_paragraph(doc,f"We modeled log ASRs using {trend_method}. Models allowed zero to two change points, required at least four annual observations per segment, used {meta['permutations']:,} residual permutations for sequential model selection, and reported segment annual percentage changes and overall average annual percentage changes (AAPCs) with model-based 95% confidence intervals. Primary models were homoscedastic because GBD UIs are not independent sampling standard errors. Weighting by log-scale standard errors approximated from marginal UIs was a sensitivity analysis.")
    add_paragraph(doc,"We tested trajectory parallelism for China versus the United States within outcome and sex, and for females versus males within outcome and country. Benjamini-Hochberg correction was applied separately to the six country and six sex comparisons. Native GBD UIs were reported for original estimates. Derived ratios, differences, and changes remained point estimates because posterior draws and cross-estimate correlations were unavailable.")
    add_heading(doc,"Demographic decomposition",2)
    add_paragraph(doc,"For each outcome, country, and sex, we decomposed the change in reconstructed adult counts into total adult population growth, changing age composition, and changing age-specific rates. We averaged marginal contributions over all six possible factor-replacement orders, which is equivalent to a Shapley decomposition. Components were required to sum to total reconstructed change within numerical tolerance. Primary estimates used 1990-2023; 2000-2023, 2010-2023, annual-chain, and five-year-chain analyses assessed endpoint and path sensitivity.")
    add_heading(doc,"Secondary age-period-cohort analysis",2)
    add_paragraph(doc,"Secondary incidence analysis used 11 equal five-year age groups from 15-19 through 65-69 years and six equal periods from 1994-1998 through 2019-2023. We estimated net drift, age-specific local drift, and nonlinear age, period, and cohort curvature contrasts. The design explicitly separated estimable curvature from the unidentified linear age-period-cohort dependency. We did not interpret period or cohort patterns causally. Excluding 2019-2023 assessed sensitivity to the pandemic-era endpoint.")
    add_heading(doc,"Reporting and reproducibility",2)
    add_paragraph(doc,"The analysis was conducted in Python with a deterministic random seed. Every table and figure is generated from saved machine-readable inputs. Reporting follows the Guidelines for Accurate and Transparent Health Estimates Reporting (GATHER). The complete code, provenance table, data dictionary, and submission-gate metadata accompany the manuscript.")

    add_heading(doc,"Results")
    add_heading(doc,"Data completeness and outcome selection",2)
    audit=tables["data_audit"]; ident=tables["yld_daly_identity"].iloc[0]; recon=tables["population_reconstruction"]
    add_paragraph(doc,f"All 12 primary country-sex-outcome panels contained 34 annual all-age counts and 34 annual ASRs. Each decomposition panel contained 12 adult age groups. No invalid UI ordering or nonpositive point estimates were identified. DALYs and YLDs were numerically identical across {int(ident.matched_cells):,} matched cells (maximum relative point-estimate difference {ident.max_relative_difference_val:.2e}); YLDs were therefore omitted as duplicate outcomes. The 99th percentile absolute population-rate reconstruction discrepancy was {recon.relative_error_pct.abs().quantile(.99):.3g}%.")
    add_heading(doc,"Burden levels and temporal changes",2)
    add_paragraph(doc,"All-age incidence, prevalence, and DALY counts increased between 1990 and 2023 in both countries and both sexes. In contrast, Chinese ASRs were nearly stable, whereas U.S. ASRs declined across the three outcomes. Males generally had higher rates than females. Thus, count growth and standardized-rate trends conveyed different dimensions of population burden.")
    table1=compact_endpoint_table(endpoints,segmented)
    add_paragraph(doc,"Table 1. Endpoint burden and modeled trends",bold_lead="Table 1.")
    add_dataframe_table(doc,table1,[("Location","Location"),("Sex","Sex"),("Outcome","Outcome"),("Count 2023 (95% UI)","Count 2023 (95% UI)"),("Count change, %","Count change, %"),("ASR 1990","ASR 1990"),("ASR 2023 (95% UI)","ASR 2023 (95% UI)"),("AAPC, %/year","AAPC, %/year")],font_size=7.2)
    add_figure(doc,analysis/"figures"/"main"/"figure_1_asr_trends.png","Figure 1. Sex-specific schizophrenia age-standardized rates in China and the United States, 1990-2023. Shading shows native 95% GBD uncertainty intervals.","Three panels showing incidence, prevalence, and DALY age-standardized rate trends by country and sex.")
    add_heading(doc,"Segmented trends and formal comparisons",2)
    rejected=pairs[pairs.q_value_bh_within_family<.05]
    add_paragraph(doc,f"Permutation-selected models identified nonlinear rate trajectories in several panels. False-discovery-rate-adjusted parallelism tests rejected parallel trends in {len(rejected)} of 12 comparisons. These tests quantify trajectory differences and do not identify causal explanations for them.")
    pair_display=pairs.copy()
    pair_display["group_a"]=pair_display.group_a.str.replace("United States of America","United States",regex=False)
    pair_display["group_b"]=pair_display.group_b.str.replace("United States of America","United States",regex=False)
    pair_display["stratum"]=pair_display.stratum.str.replace("United States of America","United States",regex=False)
    pair_display["Comparison"]=pair_display.group_a+" vs "+pair_display.group_b
    pair_display["p"]=pair_display.p_value.map(lambda x:f"{x:.3g}"); pair_display["q"]=pair_display.q_value_bh_within_family.map(lambda x:f"{x:.3g}")
    pair_display["Nonparallel"]=pair_display.q_value_bh_within_family.map(lambda x:"Yes" if x<.05 else "No")
    add_paragraph(doc,"Table 2. Prespecified trajectory parallelism tests",bold_lead="Table 2.")
    add_dataframe_table(doc,pair_display,[("comparison_family","Family"),("stratum","Stratum"),("measure_name","Outcome"),("Comparison","Comparison"),("p","P value"),("q","BH q value"),("Nonparallel","Nonparallel at q<0.05")],font_size=8)
    add_figure(doc,analysis/"figures"/"main"/"figure_2_segmented_trends.png","Figure 2. Observed and fitted segmented age-standardized rate trajectories. The provisional curves use the independent permutation implementation and are not NCI Joinpoint output.","Six panels showing observed and fitted sex-specific segmented rate trends.")
    add_heading(doc,"Age-specific patterns",2)
    add_paragraph(doc,"The age profiles differed by outcome and sex, but male rates were generally higher. Incidence peaked in younger adulthood, whereas prevalence and DALY rates remained substantial across middle adulthood. The coarse 70+ terminal category limited interpretation of late-life heterogeneity.")
    add_figure(doc,analysis/"figures"/"main"/"figure_3_age_patterns.png","Figure 3. Age-specific incidence, prevalence, and DALY rates in 2023 by country and sex. Shading shows native 95% GBD uncertainty intervals.","Six panels showing age-specific schizophrenia rates in 2023 by country and sex.")
    add_heading(doc,"Drivers of changing adult counts",2)
    primary=decomp[(decomp.start_year==1990)&(decomp.end_year==2023)].copy()
    primary["Location"]=primary.location_name.str.replace("United States of America","United States",regex=False)
    for c in ("population_growth","population_aging","rate_change","total_change"): primary[c]=primary[c].map(lambda x:fmt(x,0))
    add_paragraph(doc,"Table 3. Shapley decomposition of adult count changes, 1990-2023",bold_lead="Table 3.")
    add_dataframe_table(doc,primary,[("Location","Location"),("sex_name","Sex"),("measure_name","Outcome"),("population_growth","Population growth"),("population_aging","Population ageing"),("rate_change","Rate change"),("total_change","Total change")],font_size=8)
    add_paragraph(doc,"Population growth was the dominant positive contribution in most panels. Population ageing and age-specific rate change varied in direction and magnitude by country, sex, and outcome. Negative components represent countervailing forces; component percentages can therefore exceed 100% and were not treated as compositional shares.")
    add_figure(doc,analysis/"figures"/"main"/"figure_4_decomposition.png","Figure 4. Shapley decomposition of changes in reconstructed adult schizophrenia counts, 1990-2023. Components are deterministic attributions based on posterior mean rates.","Six panels decomposing adult count changes into population growth, ageing, and rate change.",width=5.55)
    add_heading(doc,"Secondary age-period-cohort summaries",2)
    add_paragraph(doc,"Incidence net and local drifts were directionally consistent with the primary rate-trend results. Nonlinear age, period, and cohort curvature contrasts are presented in the supplement. Their interpretation remained descriptive because the exact linear dependency among age, period, and cohort prevents unique causal separation.")

    add_heading(doc,"Discussion")
    add_paragraph(doc,"This comparative study found a shared increase in the absolute number of people affected by schizophrenia-related outcomes but divergent standardized-rate trajectories. Chinese age-standardized rates were broadly stable, whereas U.S. rates declined, particularly among females. The apparent contradiction between rising counts and flat or falling ASRs was largely explained by population growth and changing population structure.")
    add_paragraph(doc,"The distinction has direct public-health relevance. Counts approximate the scale of service demand, while standardized rates better describe changes in population risk after controlling for age structure. Service planning based only on standardized rates may underestimate future capacity needs; interpreting count growth as evidence of worsening individual risk would be equally misleading.")
    add_paragraph(doc,"Sex-specific analyses showed persistent male excess for several outcomes, consistent with established differences in age at onset, course, and disability. However, modeled GBD estimates cannot determine whether country or sex contrasts arise from true incidence, diagnostic recognition, data availability, care access, remission, excess mortality, or modeling assumptions. Health-system differences are therefore contextual hypotheses rather than tested mechanisms.")
    add_heading(doc,"Strengths and limitations",2)
    add_paragraph(doc,"Strengths include a prespecified comparative question, consistent GBD 2023 outcome definitions, explicit separation of native UIs from model CIs, formal multiplicity-controlled trajectory comparisons, exact Shapley decomposition, and restriction of APC interpretation to estimable functions. The pipeline records exclusions and prevents a provisional population source from being silently treated as official.")
    add_paragraph(doc,"Several limitations are important. First, GBD values are modeled estimates rather than direct observations and may share smoothing assumptions across years and countries. Second, posterior draws were unavailable, so exact uncertainty could not be propagated to changes, ratios, decomposition components, or comparative tests. Third, cross-year and cross-stratum correlations were unknown. Fourth, the 70+ age category obscured late-life patterns. Fifth, APC estimates remain sensitive to grouping and cannot identify independent causal age, period, and cohort effects. Sixth, ecological country contrasts cannot support causal attribution to policy or health systems. Finally, the provisional build requires the authenticated official GBD 2023 population export and user-registered NCI output before submission.")
    add_heading(doc,"Conclusions")
    add_paragraph(doc,"Between 1990 and 2023, schizophrenia incidence, prevalence, and DALY counts increased in China and the United States even as standardized-rate trajectories diverged. Population growth drove much of the increase in absolute burden, while age-specific rate change distinguished the two countries. Public-health planning should jointly consider service-volume counts, standardized epidemiologic trends, age structure, and sex-specific needs.")

    add_heading(doc,"Declarations")
    add_heading(doc,"Ethics approval and consent to participate",2); add_paragraph(doc,"Not applicable. This study used aggregated, non-identifiable modeled estimates available through the Institute for Health Metrics and Evaluation.")
    add_heading(doc,"Consent for publication",2); add_paragraph(doc,"Not applicable.")
    add_heading(doc,"Availability of data and materials",2); add_paragraph(doc,"GBD estimates are available through the IHME GBD Results Tool subject to IHME terms. Analytic code, derived tables, provenance metadata, and exact model settings accompany this submission. The repository does not redistribute data beyond applicable IHME terms.")
    add_heading(doc,"Competing interests",2); add_paragraph(doc,"The authors declare no competing interests.")
    add_heading(doc,"Funding",2); add_paragraph(doc,"No study-specific funding is declared in this blinded draft.")
    add_heading(doc,"Authors' contributions",2); add_paragraph(doc,"Contributor roles will be reported using the CRediT taxonomy in the unblinded submission file.")

    add_heading(doc,"References")
    refs=[
        "1. GBD 2021 Mental Disorders Collaborators. Global, regional, and national burden of 12 mental disorders in 204 countries and territories, 1990-2019. Lancet Psychiatry. 2022;9:137-150.",
        "2. Charlson FJ, Ferrari AJ, Santomauro DF, et al. Global epidemiology and burden of schizophrenia: findings from the Global Burden of Disease Study 2016. Schizophr Bull. 2018;44:1195-1203.",
        "3. Solmi M, Seitidis G, Mavridis D, et al. Incidence, prevalence, and global burden of schizophrenia: data, with critical appraisal, from GBD 2019. Mol Psychiatry. 2023;28:5319-5327.",
        "4. Global Burden of Disease Collaborative Network. Global Burden of Disease Study 2023 results. Institute for Health Metrics and Evaluation; 2026. https://vizhub.healthdata.org/gbd-results/.",
        "5. Stevens GA, Alkema L, Black RE, et al. Guidelines for Accurate and Transparent Health Estimates Reporting: the GATHER statement. Lancet. 2016;388:e19-e23.",
        "6. Kim HJ, Fay MP, Feuer EJ, Midthune DN. Permutation tests for joinpoint regression with applications to cancer rates. Stat Med. 2000;19:335-351; correction 2001;20:655.",
        "7. National Cancer Institute. Joinpoint Regression Program, version 6.0.1. Surveillance Research Program; 2026. https://surveillance.cancer.gov/joinpoint/.",
        "8. Das Gupta P. Standardization and decomposition of rates: a user's manual. US Bureau of the Census; 1993.",
        "9. Holford TR. The estimation of age, period and cohort effects for vital rates. Biometrics. 1983;39:311-324.",
        "10. Clayton D, Schifflers E. Models for temporal variation in cancer rates. II: age-period-cohort models. Stat Med. 1987;6:469-481.",
        "11. Luo L. Assessing validity and application scope of the intrinsic estimator approach to the age-period-cohort problem. Demography. 2013;50:1945-1967.",
        "12. Rutherford MJ, Lambert PC, Thompson JR. Age-period-cohort modeling. Stata J. 2010;10:606-627.",
    ]
    for ref in refs:
        p_ref=doc.add_paragraph(ref,style=None)
        p_ref.paragraph_format.first_line_indent=Inches(-.2); p_ref.paragraph_format.left_indent=Inches(.2)
        p_ref.paragraph_format.space_after=Pt(3)
        for run in p_ref.runs: run.font.size=Pt(9.3)

    out=out_dir/"manuscript_BMC_Public_Health.docx"; doc.save(out); return out


def build_supplement(analysis: Path, out_dir: Path, meta: dict, tables: dict[str,pd.DataFrame]):
    doc=setup_document("Supplementary material","China-US schizophrenia burden study, 1990-2023",not meta["submission_ready"])
    add_heading(doc,"Supplementary methods and audit trail")
    add_paragraph(doc,"This supplement contains the full audit, sensitivity analyses, secondary APC outputs, and external submission gates. Native GBD uncertainty intervals are never relabeled as confidence intervals for derived quantities.")
    add_heading(doc,"Table S1. Data provenance",2)
    provenance=tables["provenance"].copy()
    provenance["file"]=provenance["file"].map(lambda x:Path(str(x)).name if ("/" in str(x) or "\\" in str(x)) else str(x))
    provenance["status"]=provenance["status"].replace({"derived_proxy_NOT_OFFICIAL":"PROVISIONAL proxy"})
    add_dataframe_table(doc,provenance,[(c,c.replace("_"," ").title()) for c in provenance.columns],font_size=7.5)
    doc.add_page_break()
    add_heading(doc,"Table S2. Completeness and validity audit",2)
    add_dataframe_table(doc,tables["data_audit"],[(c,c.replace("_"," ").title()) for c in tables["data_audit"].columns],font_size=7.2)
    add_heading(doc,"Table S3. YLD-DALY identity audit",2)
    ident=tables["yld_daly_identity"].copy()
    for c in ident.columns:
        if c.startswith("max_"): ident[c]=ident[c].map(lambda v:f"{float(v):.3e}")
    add_dataframe_table(doc,ident,[(c,c.replace("_"," ").title()) for c in ident.columns],font_size=7.5)
    doc.add_page_break()
    add_heading(doc,"Table S4. Weighted trend sensitivity",2)
    x=tables["ui_weighted_sensitivity"].copy()
    for c in ("primary_aapc","ui_weighted_fixed_knot_aapc","difference"): x[c]=x[c].map(lambda v:fmt(v,4))
    add_dataframe_table(doc,x,[(c,c.replace("_"," ").title()) for c in x.columns],font_size=7.5)
    add_heading(doc,"Table S5. Alternative-window decomposition",2)
    d=tables["decomposition"][tables["decomposition"].start_year.isin([2000,2010])].copy()
    d["location_name"]=d.location_name.str.replace("United States of America","United States",regex=False)
    for c in (*("population_growth","population_aging","rate_change","total_change"),):
        if c in d: d[c]=d[c].map(lambda v:fmt(v,1))
    d["closure_error"]=d.closure_error.map(lambda v:"<1e-8" if abs(float(v))<1e-8 else f"{float(v):.2e}")
    add_dataframe_table(doc,d,[(c,c.replace("_"," ").title()) for c in ["location_name","sex_name","measure_name","start_year","end_year","population_growth","population_aging","rate_change","total_change","closure_error"]],font_size=7.5)
    add_heading(doc,"Secondary APC analysis",1)
    add_paragraph(doc,"The APC analysis is intentionally restricted to incidence and estimable summaries. Curvature relative risks describe nonlinear departures after removing intercept and linear trend. They are not separately identified causal effects.")
    apc=tables["apc_summary"].copy(); apc["location_name"]=apc.location_name.str.replace("United States of America","United States",regex=False)
    for c in ("net_drift","net_drift_lower_model_ci","net_drift_upper_model_ci"): apc[c]=apc[c].map(lambda v:f"{float(v):.4f}")
    add_dataframe_table(doc,apc,[(c,c.replace("_"," ").title()) for c in apc.columns],font_size=7.5)
    add_figure(doc,analysis/"figures"/"supplement"/"figure_s1_counts.png","Figure S1. All-age schizophrenia outcome counts by country and sex, 1990-2023.","Three panels showing all-age incidence, prevalence, and DALY counts.")
    add_figure(doc,analysis/"figures"/"supplement"/"figure_s2_apc_incidence.png","Figure S2. Secondary incidence APC estimable summaries. Relative-risk curves are nonlinear curvature contrasts, not independent causal age, period, or cohort effects.","Four panels showing incidence net drift and nonlinear age, period, and cohort curvature summaries.")
    add_heading(doc,"External submission gates")
    add_bullets(doc,[
        f"Population status: {meta['population_status']}",
        f"Trend status: {meta['trend_status']}",
        "Official GBD 2023 population must replace the proxy before submission.",
        "Official NCI Joinpoint output must be generated by a registered end user and imported before the method is described as NCI Joinpoint.",
        "Posterior draws remain unavailable; derived uncertainty must remain explicitly limited.",
    ])
    out=out_dir/"supplementary_material.docx"; doc.save(out); return out


def build_methods_appendix(out_dir: Path, meta: dict):
    doc=setup_document("Statistical methods appendix","China-US schizophrenia burden study",False,show_blinded=False)
    if not meta["submission_ready"]:
        note=doc.add_paragraph()
        run=note.add_run("Provisional methods appendix: official GBD 2023 population and registered NCI Joinpoint outputs remain required before submission.")
        run.bold=True; run.font.color.rgb=RGBColor(156,87,0)
    add_heading(doc,"Analysis estimands")
    add_bullets(doc,["All-age numbers describe service-volume burden.","Age-standardized rates describe temporal and cross-country rate patterns after standardization.","Derived ratios and changes are point estimands; no draw-based uncertainty is asserted.","Decomposition components are deterministic functions of posterior mean rates and populations."])
    add_heading(doc,"Segmented log-linear trend model")
    add_paragraph(doc,"For year t and rate r(t), log r(t) = beta0 + beta1(t-t0) + sum[j] delta_j max(0,t-tau_j) + error(t). Candidate models contain zero, one, or two knots and require four observations between boundaries. Sequential residual-permutation tests compare zero versus one and one versus two knots. Segment APC is 100[exp(slope)-1]; AAPC is the time-weighted mean log slope transformed to a percentage.")
    add_paragraph(doc,f"The provisional build uses {meta['permutations']:,} permutations and random seed {meta['seed']}. It is an independent implementation, not NCI Joinpoint. Official NCI inputs and settings are exported separately.")
    add_heading(doc,"Parallelism tests")
    add_paragraph(doc,"For each prespecified pair, a common spline basis used the union of selected knots. The restricted model allowed a group intercept but common slopes; the full model added group-by-slope interactions. An F test compared residual sums of squares. Benjamini-Hochberg correction was applied within the country-comparison and sex-comparison families.")
    add_heading(doc,"Shapley decomposition")
    add_paragraph(doc,"For total adult population N, age-share vector S, and age-specific rate vector R, reconstructed burden is B=N sum_a(S_a R_a). Each factor was replaced from baseline to endpoint in all 3! orders. A factor's contribution is the average marginal change across orders. The three contributions therefore sum exactly to B_end-B_start apart from floating-point tolerance.")
    add_heading(doc,"APC estimable functions")
    add_paragraph(doc,"The exact relation cohort=period-age prevents simultaneous identification of unrestricted linear age, period, and cohort effects. The model therefore contains an intercept, identifiable longitudinal age slope, net drift, and nonlinear age, period, and cohort bases constrained to be orthogonal to intercept and linear trend. Net and local drift are estimated separately from annual log-rate regressions. Reported relative risks are curvature contrasts to central reference categories.")
    doc.add_page_break()
    add_heading(doc,"Uncertainty taxonomy")
    add_dataframe_table(doc,pd.DataFrame([
        {"Quantity":"Original GBD estimate","Interval":"Native 95% UI","Interpretation":"2.5th-97.5th percentiles of GBD draws"},
        {"Quantity":"Segment APC/AAPC","Interval":"Model 95% CI","Interpretation":"Conditional regression uncertainty; not GBD posterior uncertainty"},
        {"Quantity":"Ratios and changes","Interval":"None","Interpretation":"Point estimate; interval dominance is conservative only"},
        {"Quantity":"Decomposition component","Interval":"None","Interpretation":"Deterministic attribution using posterior means"},
    ]),[("Quantity","Quantity"),("Interval","Reported uncertainty"),("Interpretation","Interpretation")],font_size=8.5)
    out=out_dir/"statistical_methods_appendix.docx"; doc.save(out); return out


def build_gather(out_dir: Path, provisional: bool):
    doc=setup_document("GATHER reporting checklist","Guidelines for Accurate and Transparent Health Estimates Reporting",provisional)
    rows=[
        (1,"Define indicators, populations, and time periods","Methods: Outcomes; Study design","Complete"),
        (2,"List funding sources and funder roles","Declarations: Funding","Complete for blinded draft"),
        (3,"Describe source data access","Methods: Data source; Availability statement","Complete"),
        (4,"Identify all data sources","Supplementary Table S1","Complete"),
        (5,"Describe data-source inclusion/exclusion","Methods: Prespecified exclusions","Complete"),
        (6,"Report data-source characteristics","Supplementary Table S1","Complete"),
        (7,"Provide source data in machine-readable form","Provenance CSV and source archives","Complete subject to IHME terms"),
        (8,"Describe data processing","Methods; Statistical appendix","Complete"),
        (9,"Describe model and parameter selection","Methods: Trend/APC; Statistical appendix","Complete"),
        (10,"Describe covariates","No covariates used","Not applicable"),
        (11,"Describe uncertainty analysis","Methods: Comparisons; Statistical appendix","Complete with limitation"),
        (12,"Describe analytical or statistical software","Methods: Reproducibility; README","Complete"),
        (13,"Provide source code","publication_study source files","Complete"),
        (14,"Provide results in machine-readable form","CSV and XLSX outputs","Complete"),
        (15,"Report quantitative estimates with uncertainty","Results; Table 1; Figures 1 and 3","Complete for native estimates"),
        (16,"Make comparison groups explicit","Methods and all tables","Complete"),
        (17,"Discuss limitations of data and methods","Discussion: Strengths and limitations","Complete"),
        (18,"State interpretation and implications","Discussion and Conclusions","Complete"),
    ]
    add_dataframe_table(doc,pd.DataFrame(rows,columns=["Item","Requirement","Location","Status"]),[("Item","Item"),("Requirement","Requirement"),("Location","Manuscript location"),("Status","Status")],font_size=8.5)
    if provisional: add_paragraph(doc,"Checklist completion does not override the two external submission gates identified in build_metadata.json.")
    out=out_dir/"GATHER_checklist.docx"; doc.save(out); return out


def load_tables(analysis: Path) -> dict[str,pd.DataFrame]:
    return {p.stem:pd.read_csv(p) for p in (analysis/"tables").glob("*.csv")}


def main():
    p=argparse.ArgumentParser(); p.add_argument("--analysis-dir",type=Path,required=True); p.add_argument("--output-dir",type=Path)
    args=p.parse_args(); analysis=args.analysis_dir; out=args.output_dir or analysis/"documents"; out.mkdir(parents=True,exist_ok=True)
    meta=json.loads((analysis/"build_metadata.json").read_text(encoding="utf-8")); tables=load_tables(analysis)
    paths=[build_manuscript(analysis,out,meta,tables),build_supplement(analysis,out,meta,tables),build_methods_appendix(out,meta),build_gather(out,not meta["submission_ready"])]
    (out/"document_manifest.json").write_text(json.dumps([str(x) for x in paths],indent=2),encoding="utf-8")
    print("\n".join(map(str,paths)))


if __name__=="__main__": main()
