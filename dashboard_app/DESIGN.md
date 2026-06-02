---
name: BRCA Navigator
description: Multi-omics breast cancer research dashboard with guided data exploration
colors:
  primary: "#7A70BA"
  primary-light: "#9B92CC"
  primary-dark: "#41386B"
  primary-bg: "#EDEBF5"
  semantic-protective: "#6B8F5E"
  semantic-protective-bg: "#E2EBD9"
  semantic-risk: "#B85450"
  semantic-risk-bg: "#F0DCDA"
  semantic-caution: "#B8963E"
  semantic-caution-bg: "#F3ECDA"
  neutral-paper: "#EBEED5"
  neutral-surface: "#E2E5D0"
  neutral-surface-raised: "#F5F7E8"
  neutral-border: "#C5C8B8"
  ink-primary: "#2D2A4A"
  ink-secondary: "#4A4670"
  ink-tertiary: "#6B6890"
  ink-muted: "#8E8BAE"
  sidebar: "#41386B"
  header: "#41386B"
  subtype-luminal-a: "#7A70BA"
  subtype-luminal-b: "#B8963E"
  subtype-basal: "#B85450"
  subtype-her2: "#9B6EB0"
  subtype-normal: "#6B8F5E"
typography:
  display:
    fontFamily: "Figtree, system-ui, sans-serif"
    fontSize: "2rem"
    fontWeight: 800
    lineHeight: 1.2
  headline:
    fontFamily: "Figtree, system-ui, sans-serif"
    fontSize: "1.5rem"
    fontWeight: 700
    lineHeight: 1.3
  title:
    fontFamily: "Figtree, system-ui, sans-serif"
    fontSize: "1.1rem"
    fontWeight: 600
    lineHeight: 1.4
  body:
    fontFamily: "Figtree, system-ui, sans-serif"
    fontSize: "0.95rem"
    fontWeight: 400
    lineHeight: 1.7
  label:
    fontFamily: "Figtree, system-ui, sans-serif"
    fontSize: "0.85rem"
    fontWeight: 500
    lineHeight: 1.4
    letterSpacing: "0.01em"
  caption:
    fontFamily: "Figtree, system-ui, sans-serif"
    fontSize: "0.78rem"
    fontWeight: 500
    lineHeight: 1.4
    letterSpacing: "0.03em"
rounded:
  sm: "6px"
  md: "10px"
  lg: "16px"
  full: "9999px"
spacing:
  xs: "4px"
  sm: "8px"
  md: "16px"
  lg: "24px"
  xl: "32px"
  xxl: "48px"
components:
  factor-card:
    backgroundColor: "{colors.neutral-surface}"
    textColor: "{colors.ink-primary}"
    rounded: "{rounded.md}"
    padding: "20px"
  factor-card-protective:
    backgroundColor: "{colors.semantic-protective-bg}"
    textColor: "{colors.ink-primary}"
    rounded: "{rounded.md}"
    padding: "20px"
  factor-card-risk:
    backgroundColor: "{colors.semantic-risk-bg}"
    textColor: "{colors.ink-primary}"
    rounded: "{rounded.md}"
    padding: "20px"
  simulator-card:
    backgroundColor: "{colors.neutral-surface-raised}"
    textColor: "{colors.ink-primary}"
    rounded: "{rounded.lg}"
    padding: "30px"
  simulator-result:
    backgroundColor: "{colors.primary-dark}"
    textColor: "#ffffff"
    rounded: "{rounded.lg}"
    padding: "30px"
  description-box:
    backgroundColor: "{colors.neutral-surface}"
    textColor: "{colors.ink-secondary}"
    rounded: "{rounded.md}"
    padding: "16px 20px"
  plot-container:
    backgroundColor: "{colors.neutral-surface-raised}"
    textColor: "{colors.ink-primary}"
    rounded: "{rounded.sm}"
    padding: "0"
  nav-sidebar:
    backgroundColor: "{colors.sidebar}"
    textColor: "#B1B4C8"
    rounded: "0"
    padding: "0"
  nav-sidebar-active:
    backgroundColor: "{colors.primary}"
    textColor: "#ffffff"
    rounded: "0"
    padding: "0"
---

# Design System: BRCA Navigator

## 1. Overview

**Creative North Star: "The Guided Atlas"**

This is a research tool that guides users through complex multi-omics data with the quiet authority of a well-organized atlas. Every screen serves a purpose: explore factors, discover gene associations, simulate survival outcomes. The interface recedes so the data can speak.

The system draws from a restrained palette of Delft Blue, Amethyst, and soft Beige — a calm, academic color scheme that avoids both the sterile gray of typical dashboards and the overstimulation of data-viz-heavy tools.

**Key Characteristics:**
- Single-family typography (Figtree) with weight contrast for hierarchy
- Calm, academic palette: Amethyst accent on warm beige surfaces with Delft Blue navigation
- Gently lifted depth: subtle ambient shadows, never loud
- Single theme with consistent day-one clarity
- Data-forward layout: plots and tables dominate, controls are compact

## 2. Colors

The palette is built around a calm, academic triad: Delft Blue navigation, Amethyst primary accent, and Beige content surfaces. Semantic colors are reserved exclusively for data meaning (protective, risk, caution).

### Primary
- **Amethyst** (#7A70BA): Primary accent for active states, buttons, sidebar highlights, links, and selection indicators. Used sparingly.
- **Amethyst Light** (#9B92CC): Lighter variant for hover states and subtle background tints.
- **Delft Blue** (#41386B): Header and sidebar background. The structural anchor of the palette.
- **Amethyst Ghost** (#EDEBF5): Near-white tint for background washes on cards and containers.

### Semantic (Data-Driven)
- **Protective Olivine** (#6B8F5E): Used exclusively for protective factors (HR < 1), positive outcomes, survival advantages. Background: #E2EBD9.
- **Risk Muted Red** (#B85450): Used exclusively for risk factors (HR > 1), danger signals, aggressive biology. Background: #F0DCDA.
- **Caution Warm Amber** (#B8963E): Used for warnings, borderline significance, age-related factors. Background: #F3ECDA.

### Neutral
- **Paper** (#EBEED5): Default content background. Warm beige, soft on the eyes.
- **Surface** (#E2E5D0): Elevated cards, secondary containers, slightly deeper beige.
- **Surface Raised** (#F5F7E8): Simulator cards, prominent panels. Lightest beige.
- **Border** (#C5C8B8): Standard borders on cards, inputs, dividers. Warm gray.
- **Ink Primary** (#2D2A4A): Headings, primary text. Deep navy-purple for maximum contrast on beige surfaces.
- **Ink Secondary** (#4A4670): Body text, descriptions. Comfortable reading density.
- **Ink Tertiary** (#6B6890): Labels, captions, metadata.
- **Ink Muted** (#8E8BAE): Disabled states, placeholder text, lowest-emphasis text.

### Subtype Colors
Used exclusively for PAM50 molecular subtype encoding across all plots and legends.
- **Luminal A** (#7A70BA): ~46% of patients. Best prognosis. Amethyst tint.
- **Luminal B** (#B8963E): ~18%. Intermediate prognosis. Warm amber.
- **Basal-like** (#B85450): ~15%. Most aggressive. Muted red.
- **HER2-enriched** (#9B6EB0): ~14%. Aggressive but treatable. Orchid.
- **Normal-like** (#6B8F5E): ~8%. Variable prognosis. Olivine.

### Named Rules

**The Amethyst Reserve Rule.** The primary Amethyst accent appears on no more than 10% of any given screen. Its rarity signals importance: active nav item, current selection, primary action. Everything else is neutral or semantic.

**The Semantic-Only Rule.** Green, red, and amber exist solely to encode data meaning (protective, risk, caution). They are never used for decoration, borders on empty containers, or background fills on non-data elements.

## 3. Typography

**Font:** Figtree (with system-ui, sans-serif fallback)

**Character:** A geometric sans-serif with warm terminals and generous x-height. Figtree reads as modern and approachable without being casual. The wide weight range (300-800) creates clear hierarchy through weight contrast rather than size exaggeration.

### Hierarchy
- **Display** (800, 2rem/1.2): Page-level headings in the scrollytelling narrative. Used sparingly, one per section.
- **Headline** (700, 1.5rem/1.3): Section headings within dashboard tabs. Primary information anchor.
- **Title** (600, 1.1rem/1.4): Card titles, subsection headers, stat labels. Workhorse heading.
- **Body** (400, 0.95rem/1.7): Prose, descriptions, explanatory text. Capped at 65-75ch for readability.
- **Label** (500, 0.85rem/1.4, +0.01em tracking): Form labels, nav items, buttons. Functional text.
- **Caption** (500, 0.78rem/1.4, +0.03em tracking): Table headers, metadata, fine print. Least emphasis.

### Named Rules

**The Weight Contrast Rule.** Hierarchy is established through weight (400 vs 600 vs 800), not through size alone. This keeps the type scale tight (1.125-1.2 ratio between steps) and avoids display-size sprawl.

**The Single Family Rule.** Figtree is the only typeface. The weight range is wide enough to carry all roles.

## 4. Elevation

The system uses gentle ambient shadows to lift cards and panels above the background surface. Shadows are structural (conveying hierarchy) rather than decorative.

### Shadow Vocabulary
- **Resting** (`box-shadow: 0 2px 8px rgba(0,0,0,0.06)`): Factor cards, description boxes, subtle containers at rest.
- **Elevated** (`box-shadow: 0 4px 20px rgba(0,0,0,0.08)`): Simulator cards, prominent panels, modal-like surfaces.
- **Plot** (`box-shadow: 0 1px 4px rgba(0,0,0,0.08)`): Plot containers, low-profile elevated surfaces.
- **Hover** (`box-shadow: 0 8px 24px rgba(0,0,0,0.12)`): Interactive elements on hover state.

### Named Rules

**The Flat-By-Default Rule.** Content areas (prose, tables, plain backgrounds) have no shadow. Shadows are reserved for cards, panels, and containers that sit above the surface layer.

## 5. Components

### Factor Cards
- **Character:** Information-dense summary tiles for each MOFA factor. Three variants: neutral (surface), protective (olivine ghost), risk (muted red ghost).
- **Corner Style:** Gently curved (10px radius)
- **Background:** Neutral surface (#E2E5D0) or semantic ghost backgrounds
- **Shadow Strategy:** Resting shadow at rest, elevated on hover
- **Left Accent:** 4px solid border in factor-specific color (Amethyst, Olivine, or Muted Red)
- **Internal Padding:** 20px

### Simulator Card
- **Character:** The interactive control panel for the survival simulator. Contains sliders, action buttons, and preset profile buttons.
- **Corner Style:** Large radius (16px)
- **Background:** Surface raised (#F5F7E8)
- **Shadow Strategy:** Elevated shadow
- **Internal Padding:** 30px

### Simulator Result Panel
- **Character:** The dark panel showing predicted survival statistics (5-year survival, median, hazard ratio).
- **Corner Style:** Large radius (16px)
- **Background:** Delft Blue (#41386B)
- **Text:** White, high contrast
- **Internal Padding:** 30px

### Description Boxes
- **Character:** Page-level context boxes at the top of each tab. Contain a single paragraph explaining what the view shows.
- **Corner Style:** Gently curved (10px)
- **Background:** Surface (#E2E5D0)
- **Border:** None
- **Internal Padding:** 16px horizontal, 20px vertical

### Plot Containers
- **Character:** Wrappers around Plotly and ggplot outputs. Minimal chrome, letting the visualization speak.
- **Corner Style:** Subtle curve (6px)
- **Background:** Surface raised (#F5F7E8)
- **Border:** 1px solid #C5C8B8
- **Shadow:** Plot shadow (low-profile)

### Navigation Sidebar
- **Character:** Icon-only collapsed (60px), expands to full labels on hover (280px). Delft Blue background with Amethyst active state.
- **Background:** Delft Blue (#41386B)
- **Text:** French Gray (#B1B4C8)
- **Active State:** Amethyst (#7A70BA) background, white text
- **Hover State:** Lighter Delft Blue tint
- **Typography:** Label weight (500), 0.85rem
- **Icons:** Font Awesome 6 Free, 1.1rem

### Buttons (Action)
- **Primary:** Delft Blue (#41386B) background, white text, pill shape (50px radius), 18px 48px padding
- **Outline:** Transparent background, Delft Blue border, Delft Blue text. Hover fills with Delft Blue.
- **Success:** Protective Olivine (#6B8F5E) background, white text
- **Danger:** Risk Muted Red (#B85450) background, white text
- **Transitions:** 200ms ease-out for background and transform

### Inputs / Selectors
- **Style:** Standard Shiny inputs with Figtree font. 1px solid border (#C5C8B8), 6px radius.
- **Focus:** Amethyst border accent, subtle glow
- **Background:** Surface raised (#F5F7E8)

### Data Tables (DT)
- **Header:** Surface background (#E2E5D0), ink-primary text, bold
- **Rows:** Alternating surface/surface-raised for scanability
- **Hover:** Light Amethyst tint
- **Border:** 1px solid #C5C8B8
- **Font Size:** 0.85rem for density

### Analysis Flow Steps
- **Character:** Numbered pipeline steps with colored circle indicators.
- **Number Circle:** 40px diameter, filled with lavender (#D5D0E8), Delft Blue (#41386B) text
- **Container:** Surface-raised background, 1px border, 12px radius
- **Internal Padding:** 20px
- **Layout:** Horizontal flex with gap between number, title, badge, and content

## 6. Do's and Don'ts

### Do:
- **Do** use Amethyst (#7A70BA) as the single accent color. Its rarity signals importance.
- **Do** encode data meaning with semantic colors: olivine = protective, muted red = risk, warm amber = caution.
- **Do** use weight contrast (400 vs 600 vs 800) for type hierarchy instead of size exaggeration.
- **Do** maintain the 65-75ch line length for prose sections.
- **Do** use subtle ambient shadows (0 2px 8px rgba(0,0,0,0.06)) on cards and panels.
- **Do** keep the sidebar collapsed to icon-only (60px) by default. Expand on hover for discoverability.
- **Do** use Figtree exclusively across all surfaces.
- **Do** maintain the beige-and-deep-navy contrast for a calm, readable reading environment.

### Don't:
- **Don't** use cold gray-blue neutrals. The palette is warm beige, not cool slate.
- **Don't** use dark mode. This is a single-theme design with consistent day-one clarity.
- **Don't** overstimulate: no neon indicators, no gradient everything, no badges on badges.
- **Don't** use sci-fi aesthetics: no neon on black, no glowing borders.
- **Don't** use `border-left` greater than 1px as a colored accent stripe on cards. The 4px left accent on factor cards is the sole exception.
- **Don't** use gradient text. Use a single solid color.
- **Don't** use glassmorphism decoratively.
- **Don't** use `!important` unless absolutely necessary to override framework specificity.
- **Don't** use inline styles for repeated component patterns. Extract to CSS classes.
- **Don't** ship components with only a default state. Every interactive element needs hover, focus, active, and disabled states.
