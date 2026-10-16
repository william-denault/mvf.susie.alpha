from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    CondPageBreak,
    Frame,
    KeepTogether,
    ListFlowable,
    ListItem,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "split_variational_inference_binomial_draft.pdf"


pdfmetrics.registerFont(TTFont("DVS", r"C:\Windows\Fonts\DejaVuSerif.ttf"))
pdfmetrics.registerFont(TTFont("DVS-Bold", r"C:\Windows\Fonts\DejaVuSerif-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DVS-Italic", r"C:\Windows\Fonts\DejaVuSerif-Italic.ttf"))
pdfmetrics.registerFont(TTFont("DVMath", r"C:\Windows\Fonts\DejaVuMathTeXGyre.ttf"))
pdfmetrics.registerFontFamily(
    "DVS", normal="DVS", bold="DVS-Bold", italic="DVS-Italic"
)


PAGE_W, PAGE_H = LETTER
LEFT = RIGHT = 0.82 * inch
TOP = 0.78 * inch
BOTTOM = 0.7 * inch


def draw_page(canvas, doc):
    canvas.saveState()
    page = canvas.getPageNumber()
    canvas.setFont("DVS", 8)
    canvas.setFillColor(colors.HexColor("#5B6472"))
    canvas.drawCentredString(PAGE_W / 2, 0.36 * inch, str(page))
    canvas.restoreState()


doc = BaseDocTemplate(
    str(OUTPUT),
    pagesize=LETTER,
    leftMargin=LEFT,
    rightMargin=RIGHT,
    topMargin=TOP,
    bottomMargin=BOTTOM,
    title="A splitting approach for variational inference in non-Gaussian latent Gaussian models",
    author="Dongyue Xie, Peter Carbonetto, Matthew Stephens",
)
frame = Frame(LEFT, BOTTOM, PAGE_W - LEFT - RIGHT, PAGE_H - TOP - BOTTOM, id="body")
doc.addPageTemplates(PageTemplate(id="manuscript", frames=frame, onPage=draw_page))


styles = getSampleStyleSheet()
body = ParagraphStyle(
    "Body",
    parent=styles["BodyText"],
    fontName="DVS",
    fontSize=9.3,
    leading=13.2,
    alignment=TA_JUSTIFY,
    spaceAfter=7,
    textColor=colors.HexColor("#1F252D"),
    allowWidows=0,
    allowOrphans=0,
)
abstract_style = ParagraphStyle(
    "Abstract",
    parent=body,
    fontSize=9.1,
    leading=12.8,
    leftIndent=0.26 * inch,
    rightIndent=0.26 * inch,
)
h1 = ParagraphStyle(
    "H1",
    parent=body,
    fontName="DVS-Bold",
    fontSize=15,
    leading=18,
    alignment=TA_LEFT,
    spaceBefore=12,
    spaceAfter=7,
    keepWithNext=True,
    textColor=colors.HexColor("#203A5F"),
)
h2 = ParagraphStyle(
    "H2",
    parent=body,
    fontName="DVS-Bold",
    fontSize=11.2,
    leading=14,
    spaceBefore=9,
    spaceAfter=5,
    keepWithNext=True,
    textColor=colors.HexColor("#2C4E76"),
)
title_style = ParagraphStyle(
    "Title",
    parent=body,
    fontName="DVS-Bold",
    fontSize=21,
    leading=26,
    alignment=TA_CENTER,
    spaceAfter=20,
    textColor=colors.HexColor("#172A46"),
)
author_style = ParagraphStyle(
    "Authors",
    parent=body,
    fontSize=10.5,
    leading=14,
    alignment=TA_CENTER,
    spaceAfter=8,
)
small_center = ParagraphStyle(
    "SmallCenter",
    parent=body,
    fontSize=8.2,
    leading=11,
    alignment=TA_CENTER,
    textColor=colors.HexColor("#5B6472"),
)
eq_style = ParagraphStyle(
    "Equation",
    parent=body,
    fontName="DVMath",
    fontSize=9.6,
    leading=14.5,
    alignment=TA_CENTER,
    leftIndent=0.25 * inch,
    rightIndent=0.25 * inch,
    spaceBefore=4,
    spaceAfter=8,
    borderColor=colors.HexColor("#D9DFE8"),
    borderWidth=0.55,
    borderPadding=6,
    backColor=colors.HexColor("#F7F9FC"),
)
caption = ParagraphStyle(
    "Caption",
    parent=body,
    fontSize=8,
    leading=10.5,
    alignment=TA_LEFT,
    spaceBefore=3,
    spaceAfter=5,
)
reference = ParagraphStyle(
    "Reference",
    parent=body,
    fontSize=8,
    leading=10.8,
    leftIndent=0.18 * inch,
    firstLineIndent=-0.18 * inch,
    spaceAfter=4,
)


story = []


def P(text, style=body):
    story.append(Paragraph(text, style))


def H1(text):
    story.append(Paragraph(text, h1))


def H2(text):
    story.append(Paragraph(text, h2))


def EQ(text):
    story.append(Paragraph(text, eq_style))


def bullets(items, numbered=False):
    if numbered:
        flow = [ListItem(Paragraph(item, body), leftIndent=12) for item in items]
        story.append(
            ListFlowable(
                flow,
                bulletType="1",
                start="1",
                leftIndent=22,
                bulletFontName="DVS",
                bulletFontSize=8.7,
                spaceAfter=6,
            )
        )
        return

    rows = [[Paragraph("&#8226;", body), Paragraph(item, body)] for item in items]
    bullet_table = Table(
        rows,
        colWidths=[0.18 * inch, PAGE_W - LEFT - RIGHT - 0.18 * inch],
        hAlign="LEFT",
    )
    bullet_table.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 1),
            ]
        )
    )
    story.append(bullet_table)
    story.append(Spacer(1, 4))


# Title and abstract
story.append(Spacer(1, 0.22 * inch))
story.append(
    Table(
        [[""]],
        colWidths=[PAGE_W - LEFT - RIGHT],
        rowHeights=[4],
        style=TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#203A5F"))]),
    )
)
story.append(Spacer(1, 0.18 * inch))
P(
    "A splitting approach for variational inference in<br/>"
    "non-Gaussian latent Gaussian models",
    title_style,
)
P("Dongyue Xie&nbsp;&nbsp;&nbsp; Peter Carbonetto&nbsp;&nbsp;&nbsp; Matthew Stephens", author_style)
P("Revised working draft - 28 August 2026", small_center)
story.append(Spacer(1, 0.12 * inch))
P("<b>Abstract</b>", ParagraphStyle("AbstractHead", parent=h2, alignment=TA_CENTER, fontSize=12))
P(
    "Many analyses of high-dimensional biological data begin by transforming "
    "non-Gaussian observations and then applying a Gaussian method. Examples "
    "include log-transforming RNA sequencing counts before factor analysis and "
    "transforming accessibility proportions before functional regression. These "
    "transformations are convenient, but they discard information about the "
    "observation process and can be sensitive to arbitrary choices such as a "
    "pseudocount. We develop a modular variational framework that couples an "
    "exponential-family observation model to an existing empirical Bayes Gaussian "
    "method through a latent Gaussian splitting variable. The algorithm alternates "
    "between independent low-dimensional variational updates for the non-Gaussian "
    "observations and a Gaussian update delegated to software for normal means, "
    "regression, smoothing, matrix factorization, multivariate analysis, or "
    "structured variable selection.",
    abstract_style,
)
P(
    "We give explicit updates for Poisson and binomial observations. The binomial "
    "case is especially useful when a response aggregates binary measurements: the "
    "number of successes and trials are retained, so observations with different "
    "sampling depths contribute different information. We clarify an important "
    "limitation. With one Bernoulli observation per latent mean, the splitting "
    "variance can be weakly identified or non-identifiable; replication within a "
    "binomial observation supplies extra-binomial information. We conclude with a "
    "binomial functional regression model for single-molecule chromatin accessibility "
    "data, using wavelet-based functional fine-mapping on the latent log-odds scale.",
    abstract_style,
)
story.append(PageBreak())


H1("1. Introduction")
P(
    "Consider observations <i>y</i><sub>i</sub> whose distribution belongs to an "
    "exponential family, with latent canonical or transformed mean "
    "<i>mu</i><sub>i</sub>. Depending on the structure assigned to the joint prior "
    "of the latent means, this formulation includes empirical Bayes normal means, "
    "sparse regression, smoothing, multivariate shrinkage, and matrix factorization. "
    "The Gaussian versions of these problems are computationally mature, whereas "
    "replacing Gaussian observations by counts or proportions often requires a new "
    "inference algorithm for every structured prior."
)
EQ("y<sub>i</sub> | mu<sub>i</sub> ~ D{h(mu<sub>i</sub>)},&nbsp;&nbsp;&nbsp; mu ~ f. &nbsp;&nbsp;&nbsp; (1)")
P(
    "Our goal is to reuse Gaussian methods without first transforming the data. We "
    "introduce a latent Gaussian layer that separates the likelihood from the "
    "structured prior. One module handles the non-Gaussian likelihood; the other "
    "calls an empirical Bayes Gaussian method (EBG). This is analogous in spirit to "
    "splitting methods in optimization, where a coupled problem is replaced by two "
    "simpler problems connected through low-dimensional posterior moments."
)
P(
    "The original motivation was dimension reduction for single-cell RNA sequencing "
    "counts. A second important use arises when binary measurements are aggregated "
    "within a biological sample. In Fiber-seq, accessibility is observed on individual "
    "chromatin fibers but genetic covariates are measured at the donor level. The "
    "relevant observation is a number of accessible fibers out of a known number "
    "assayed. A proportion alone cannot distinguish one out of two from one hundred "
    "out of two hundred, and treating fibers as independent genetic samples would be "
    "pseudoreplication."
)

table_data = [
    [Paragraph("<b>Structure for b</b>", caption), Paragraph("<b>Model</b>", caption), Paragraph("<b>Gaussian engine</b>", caption)],
    ["b_i iid ~ g", "univariate mean inference", "ash, ebnm"],
    ["b = X beta, beta ~ g", "regression", "sparse EB regression"],
    ["b spatially smooth", "smoothing", "wavelet or GP EB"],
    ["b = L F'", "matrix factorization", "flash"],
    ["b_i ~ mixture N(0,U_k)", "multivariate means", "mash"],
    ["b = X B, structured B", "functional selection", "fSuSiE"],
]
t = Table(table_data, colWidths=[2.12 * inch, 1.82 * inch, 1.55 * inch], repeatRows=1)
t.setStyle(
    TableStyle(
        [
            ("FONTNAME", (0, 0), (-1, -1), "DVS"),
            ("FONTSIZE", (0, 0), (-1, -1), 7.5),
            ("LEADING", (0, 0), (-1, -1), 9.5),
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#E9EEF5")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#203A5F")),
            ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#B7BDC7")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F8FAFC")]),
            ("LEFTPADDING", (0, 0), (-1, -1), 5),
            ("RIGHTPADDING", (0, 0), (-1, -1), 5),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ]
    )
)
story.append(KeepTogether([Paragraph("Table 1. Examples of reusable Gaussian modules.", caption), t]))


H1("2. Split variational inference")
H2("2.1 The augmented model")
P("Introduce a splitting variable <i>b</i> and write")
EQ(
    "y<sub>i</sub> | mu<sub>i</sub> ~ D{h(mu<sub>i</sub>)},<br/>"
    "mu<sub>i</sub> | b<sub>i</sub>, tau<sub>i</sub><super>2</super> ~ "
    "N(b<sub>i</sub>, tau<sub>i</sub><super>2</super>),&nbsp;&nbsp;&nbsp; b ~ g. &nbsp;&nbsp;&nbsp; (2)"
)
P(
    "The splitting variance may be an algorithmic device connecting the likelihood "
    "to a Gaussian engine, or part of the scientific model representing residual "
    "variation on the latent scale. This distinction matters for identifiability. "
    "We use q(mu,b) = {product_i q(mu_i)} q(b), with each local factor "
    "q(mu_i) = N(m_i,v_i). Let bbar_i and Vb_i denote the posterior mean and "
    "variance returned by the Gaussian module."
)
EQ(
    "F(q,g;tau<super>2</super>) = sum<sub>i</sub> E log p(y<sub>i</sub>|mu<sub>i</sub>) "
    "+ sum<sub>i</sub> E log p(mu<sub>i</sub>|b<sub>i</sub>,tau<sub>i</sub><super>2</super>) "
    "+ E log{g(b)/q(b)} - sum<sub>i</sub> E log q(mu<sub>i</sub>). &nbsp;&nbsp;&nbsp; (3)"
)

story.append(CondPageBreak(1.35 * inch))
H2("2.2 Modular coordinate-ascent algorithm")
bullets(
    [
        "<b>Non-Gaussian step.</b> Update every q(mu<sub>i</sub>) using its likelihood, current prior mean bbar<sub>i</sub>, and splitting variance tau<super>2</super>.",
        "<b>Gaussian step.</b> Treat m = (m<sub>1</sub>,...,m<sub>M</sub>) as a Gaussian response with residual variance tau<super>2</super>; update q(b) and the hyperparameters of g with EBG.",
        "<b>Variance step.</b> If the splitting variance is estimated, include uncertainty from both variational blocks in its M-step.",
    ],
    numbered=True,
)
EQ(
    "tau<super>2</super> &larr; (1/M) sum<sub>i</sub> "
    "{(m<sub>i</sub>-bbar<sub>i</sub>)<super>2</super> + v<sub>i</sub> + Vb<sub>i</sub>}. &nbsp;&nbsp;&nbsp; (4)"
)
P(
    "The term Vb<sub>i</sub> is essential. Omitting it generally underestimates "
    "latent residual variation. The same update applies entrywise to matrices. "
    "Known heterogeneous splitting variances can also be used, but then the Gaussian "
    "engine must accept heterogeneous standard errors."
)

H2("2.3 Relation to a collapsed objective")
P(
    "Let f<sub>tau</sub>(mu) be the Gaussian convolution of a prior G. Its log "
    "curvature has a simple lower bound:"
)
EQ(
    "d<super>2</super> log f<sub>tau</sub>(mu) / dmu<super>2</super> = "
    "-1/tau<super>2</super> + Var(b|mu)/tau<super>4</super> "
    "&ge; -1/tau<super>2</super>. &nbsp;&nbsp;&nbsp; (5)"
)
P(
    "For q(mu)=N(m,v), this gives E log f<sub>tau</sub>(mu) "
    "at least log f<sub>tau</sub>(m) - v/(2 tau<super>2</super>). The right-hand "
    "side is exactly the result of profiling the splitting ELBO over q(b). Thus the "
    "split objective is a lower bound on the variational objective obtained after "
    "integrating out b. Equality holds when G is a point mass; with a flexible G, "
    "the gap is the price of modularity."
)


H1("3. Variational Gaussian updates")
H2("3.1 Poisson observations")
P(
    "For y<sub>i</sub>|mu<sub>i</sub> ~ Poisson{s<sub>i</sub> exp(mu<sub>i</sub>)}, "
    "the local objective for q(mu<sub>i</sub>)=N(m<sub>i</sub>,v<sub>i</sub>) is"
)
EQ(
    "F<sub>i</sub> = y<sub>i</sub>m<sub>i</sub> - s<sub>i</sub>exp(m<sub>i</sub>+v<sub>i</sub>/2) "
    "- {(m<sub>i</sub>-bbar<sub>i</sub>)<super>2</super>+v<sub>i</sub>}/(2tau<super>2</super>) "
    "+ (1/2)log v<sub>i</sub>. &nbsp;&nbsp;&nbsp; (6)"
)
P("The negative objective is strictly convex. Its stationary conditions are")
EQ(
    "0 = y<sub>i</sub> - s<sub>i</sub>exp(m<sub>i</sub>+v<sub>i</sub>/2) "
    "- (m<sub>i</sub>-bbar<sub>i</sub>)/tau<super>2</super>,<br/>"
    "1/v<sub>i</sub> = s<sub>i</sub>exp(m<sub>i</sub>+v<sub>i</sub>/2) + 1/tau<super>2</super>. &nbsp;&nbsp;&nbsp; (7)"
)
P(
    "Damped Newton updates or a one-dimensional root solver are stable. The second "
    "condition implies 0 &lt; v<sub>i</sub> &lt; tau<super>2</super>."
)

H2("3.2 Binomial observations")
P(
    "Let k<sub>i</sub>|mu<sub>i</sub> ~ Binomial{n<sub>i</sub>,logit<super>-1</super>(mu<sub>i</sub>)}, "
    "with 0 &le; k<sub>i</sub> &le; n<sub>i</sub>. The local objective is"
)
EQ(
    "F<sub>i</sub> = k<sub>i</sub>m<sub>i</sub> - n<sub>i</sub>A(m<sub>i</sub>,v<sub>i</sub>) "
    "- {(m<sub>i</sub>-bbar<sub>i</sub>)<super>2</super>+v<sub>i</sub>}/(2tau<super>2</super>) "
    "+ (1/2)log v<sub>i</sub>,<br/>"
    "A(m,v) = E<sub>Z~N(m,v)</sub> log(1+exp Z). &nbsp;&nbsp;&nbsp; (8)"
)
P("Gauss-Hermite quadrature accurately evaluates the logistic-normal expectation:")
EQ(
    "A(m,v) &asymp; (1/sqrt(pi)) sum<sub>r=1</sub><super>R</super> "
    "w<sub>r</sub> log{1+exp(m+sqrt(2v) z<sub>r</sub>)}. &nbsp;&nbsp;&nbsp; (9)"
)
P(
    "Stable code uses a softplus evaluation. Exact derivative identities are "
    "A<sub>m</sub> = E p(Z) and A<sub>v</sub> = (1/2)E p(Z){1-p(Z)}, where "
    "p(z)=logit<super>-1</super>(z). The local first-order conditions are"
)
EQ(
    "0 = k<sub>i</sub> - n<sub>i</sub>E p(Z<sub>i</sub>) "
    "- (m<sub>i</sub>-bbar<sub>i</sub>)/tau<super>2</super>,<br/>"
    "1/v<sub>i</sub> = 1/tau<super>2</super> + n<sub>i</sub>E[p(Z<sub>i</sub>){1-p(Z<sub>i</sub>)}]. &nbsp;&nbsp;&nbsp; (10)"
)
P(
    "We recommend optimization in (m<sub>i</sub>, log v<sub>i</sub>) with a damped "
    "Newton or trust-region method, using the same quadrature for the objective and "
    "derivatives. Warm starts make these independent two-parameter updates inexpensive."
)

H2("3.3 Bernoulli identifiability and replication")
P(
    "The binomial extension requires more care than the Poisson extension. Suppose "
    "x|mu ~ Bernoulli{logit<super>-1</super>(mu)} and mu ~ N(0,tau<super>2</super>). "
    "Symmetry and p(z)+p(-z)=1 imply Pr(x=1)=1/2 for every tau<super>2</super>. "
    "The Bernoulli likelihood therefore contains no information about the splitting "
    "variance in this case. Nonzero or varying prior means can supply information, "
    "but it may remain weak, and a flexible prior scale can be confounded with tau."
)
P("Replication within a binomial observation changes the variance:")
EQ(
    "Var(K) = n pbar(1-pbar) + n(n-1) Var(P),&nbsp;&nbsp;&nbsp; pbar = E(P). &nbsp;&nbsp;&nbsp; (11)"
)
P(
    "The extra-binomial term vanishes at n=1 and grows quadratically with n. Aggregated "
    "binomial data can therefore identify latent heterogeneity through repeated trials, "
    "whereas a single Bernoulli outcome cannot. Practical identification also depends "
    "on the number of independent biological units, range of trial totals, range of "
    "probabilities, and flexibility of the structured prior."
)
bullets(
    [
        "Retain trial totals and diagnose their distribution instead of analyzing proportions alone.",
        "Compare fixed and estimated splitting variances; profile the objective for flat regions or boundary optima.",
        "Use null simulations or permutations that preserve trial totals and dependence structure.",
    ]
)

H2("3.4 Alternative augmentations for binary data")
P(
    "With one Bernoulli trial per latent mean, a fixed splitting variance can still be "
    "a computational approximation but should not be interpreted as estimated "
    "overdispersion. Albert-Chib augmentation is exact for a <i>probit</i> likelihood, "
    "Pr(x=1)=Phi(b), not for the logistic sigmoid. Polya-Gamma augmentation is "
    "conditionally Gaussian for logistic binomial models. It produces observation-specific "
    "precisions, requiring a weighted Gaussian engine. Split VA trades this downstream "
    "complexity for independent local updates and a common latent residual variance."
)


H1("4. Reusable Gaussian modules")
H2("4.1 Empirical Bayes normal means")
P(
    "When b<sub>i</sub> are iid from g, the Gaussian step is the empirical Bayes normal "
    "means problem m<sub>i</sub>|b<sub>i</sub> ~ N(b<sub>i</sub>,tau<super>2</super>). "
    "Scale mixtures of normals, unimodal families, or nonparametric priors can be fitted "
    "with existing software. Their posterior means and variances feed the next split-VA "
    "iteration."
)
H2("4.2 Regression, smoothing, and multivariate means")
P(
    "For b=X beta, the Gaussian module is ordinary empirical Bayes regression. For a "
    "smooth latent function, it may be a Gaussian-process or wavelet smoother. For "
    "multivariate effects, it may fit a mixture b<sub>i</sub> ~ sum<sub>k</sub> "
    "pi<sub>k</sub>N(0,U<sub>k</sub>). The Poisson or binomial observation code remains "
    "unchanged across these choices."
)
story.append(CondPageBreak(1.35 * inch))
H2("4.3 Matrix factorization")
EQ(
    "y<sub>ij</sub>|mu<sub>ij</sub> ~ D{h(mu<sub>ij</sub>)},&nbsp;&nbsp;&nbsp; "
    "mu<sub>ij</sub>|L,F ~ N(sum<sub>r</sub> l<sub>ir</sub>f<sub>jr</sub>, tau<sub>ij</sub><super>2</super>). &nbsp;&nbsp;&nbsp; (12)"
)
P(
    "The non-Gaussian module updates each matrix entry, while the Gaussian module fits "
    "factor analysis to the matrix of latent means M=(m<sub>ij</sub>). Sparsity, "
    "non-negativity, or smoothness can be imposed through priors on loadings or factors "
    "without changing the observation module."
)

H2("4.4 Binomial functional regression and Fiber-seq")
P(
    "Let A<sub>irt</sub> be the binary accessibility call for fiber r from donor i at "
    "position t. Set K<sub>it</sub> = sum<sub>r</sub>A<sub>irt</sub>, and let "
    "N<sub>it</sub> be the number of observed fibers. A donor-level functional model is"
)
EQ(
    "K<sub>it</sub>|mu<sub>it</sub> ~ Binomial{N<sub>it</sub>,logit<super>-1</super>(mu<sub>it</sub>)},<br/>"
    "mu<sub>it</sub>|a<sub>t</sub>,B,tau<super>2</super> ~ "
    "N{a<sub>t</sub>+(XB)<sub>it</sub>,tau<super>2</super>}. &nbsp;&nbsp;&nbsp; (13)"
)
P(
    "X contains donor-level genotypes and covariates, a<sub>t</sub> is baseline "
    "accessibility, and rows of B are SNP effect curves on the log-odds scale. A "
    "functional Sum of Single Effects prior can be placed on B: each component selects "
    "one SNP, while a sparse wavelet prior regularizes its effect curve. The likelihood "
    "uses the binary sampling process, and the wavelet transform is applied to latent "
    "effect functions rather than directly to binary observations."
)
P("One binomial fSuSiE iteration would:")
bullets(
    [
        "update every q(mu<sub>it</sub>) using the binomial VGA with prior mean a<sub>t</sub> + E(XB)<sub>it</sub> and variance tau<super>2</super>;",
        "fit Gaussian fSuSiE to M=(m<sub>it</sub>), holding its residual variance fixed at tau<super>2</super>;",
        "update the baseline curve and the first two posterior moments of a<sub>t</sub>+(XB)<sub>it</sub>;",
        "update tau<super>2</super> with posterior uncertainty from the latent and fSuSiE blocks.",
    ],
    numbered=True,
)
P(
    "The raw inputs must include both K<sub>it</sub> and N<sub>it</sub>. Supplying only "
    "their ratio loses sampling precision. If fibers are more correlated than the "
    "logistic-normal layer accommodates, donor- or region-specific random effects can "
    "be added. Validation should reproduce observed coverage, baseline accessibility, "
    "within-fiber spatial dependence, and linkage disequilibrium among SNPs."
)


H1("5. Computation and diagnostics")
P(
    "The local Poisson and binomial updates are embarrassingly parallel. Their cost is "
    "linear in the number of observations times the number of quadrature nodes; the "
    "Gaussian module usually dominates running time. Recommended diagnostics include:"
)
bullets(
    [
        "monotonicity of the full split ELBO, evaluated with the same quadrature used for updates;",
        "relative changes in latent means, structured means, and the splitting variance;",
        "frequency of local solver failures and boundary values of v<sub>i</sub>;",
        "profiles of the objective over tau<super>2</super>;",
        "calibration under null simulations or covariate permutations;",
        "sensitivity to quadrature order, initialization, and fixed versus estimated tau<super>2</super>.",
    ]
)
P(
    "For binomial functional regression, observed and fitted success counts should also "
    "be compared across coverage strata and positions. A method that fits mean "
    "accessibility but fails to reproduce donor-to-donor variance is not adequately "
    "modeling the experiment, even if its fitted mean curves look plausible."
)


H1("6. Discussion")
P(
    "Split variational inference makes non-Gaussian modeling modular. The likelihood is "
    "handled by simple local variational problems, while structural learning is "
    "delegated to adaptive shrinkage, multivariate empirical Bayes, factor analysis, "
    "or functional fine-mapping. The binomial extension contributes both a useful "
    "application and a necessary warning: aggregated binary data with known totals are "
    "binomial counts whose precision depends on trial number and success probability, "
    "whereas pure Bernoulli data can provide little information about an extra latent "
    "variance."
)
P(
    "For Fiber-seq fine-mapping, the binomial formulation is scientifically more "
    "defensible than applying a homoskedastic Gaussian model to mean accessibility. It "
    "is also computationally realistic because the Gaussian fSuSiE engine can be reused "
    "once it accepts a latent response and a fixed residual variance. The main research "
    "tasks are the binomial local solver, propagation of posterior uncertainty, "
    "identifiability diagnostics, and realistic calibration experiments rather than a "
    "redesign of functional fine-mapping."
)


H1("References")
refs = [
    "Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and polychotomous response data. <i>Journal of the American Statistical Association</i>, 88:669-679.",
    "Blei, D. M., Kucukelbir, A., and McAuliffe, J. D. (2017). Variational inference: A review for statisticians. <i>Journal of the American Statistical Association</i>, 112:859-877.",
    "Boyd, S., Parikh, N., Chu, E., Peleato, B., and Eckstein, J. (2011). Distributed optimization and statistical learning via the alternating direction method of multipliers. <i>Foundations and Trends in Machine Learning</i>, 3:1-122.",
    "Denault, W. R. P. et al. (2025). fSuSiE enables fine-mapping of QTLs from genome-scale molecular profiles. <i>bioRxiv</i>, doi:10.1101/2025.08.17.670732.",
    "Jaakkola, T. S. and Jordan, M. I. (1997). A variational approach to Bayesian logistic regression models and their extensions. <i>Proceedings of the Sixth International Workshop on Artificial Intelligence and Statistics</i>, 283-294.",
    "Polson, N. G., Scott, J. G., and Windle, J. (2013). Bayesian inference for logistic models using Polya-Gamma latent variables. <i>Journal of the American Statistical Association</i>, 108:1339-1349.",
    "Stephens, M. (2017). False discovery rates: a new deal. <i>Biostatistics</i>, 18:275-294.",
    "Wang, G., Sarkar, A., Carbonetto, P., and Stephens, M. (2020). A simple new approach to variable selection in regression, with application to genetic fine mapping. <i>Journal of the Royal Statistical Society: Series B</i>, 82:1273-1300.",
]
for i, ref in enumerate(refs, 1):
    story.append(Paragraph(f"[{i}] {ref}", reference))


doc.build(story)
print(OUTPUT)
