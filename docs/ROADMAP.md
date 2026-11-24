# Fluxsheet – Roadmap und Arbeitsteilung (Claude / ChatGPT)

Stand: 2026-11-25 · Basis: `main` @ `23a8204`

Diese Datei verteilt die offenen Features auf zwei Agenten, die **gleichzeitig** arbeiten:
**Claude** (Claude Code) und **ChatGPT**. Die Pakete sind so geschnitten, dass in jeder
Welle die beiden Agenten möglichst **verschiedene Dateien** anfassen. Eine Welle ist fertig,
wenn beide Pakete in `main` gemergt sind; erst dann beginnt die nächste.

Status-Werte: `offen` · `in Arbeit` · `Review` · `gemergt`. Wer ein Paket beginnt oder
abschließt, ändert nur die Statuszeile **seines** Pakets.

**Nächster Start ist G10 in Welle 5.** Welle 5 hat Vorrang vor allen weiteren
Implementierungsarbeiten. G9 bis G12 werden in dieser Reihenfolge umgesetzt, abgenommen
und gemergt. Erst danach beginnen weitere Roadmap-Pakete oder neue Apparatefeatures.

---

## Spielregeln für paralleles Arbeiten

1. **Eigener Branch, nie direkt auf `main`.**
   - Claude: eigener Git-Worktree unter `.claude/worktrees/…`, Branch `claude/<paket-id>`.
   - ChatGPT: im Haupt-Checkout `I:\AspenPlus_Copy`, Branch `chatgpt/<paket-id>`
     (vor dem Start `git switch -c chatgpt/<paket-id> main`).
   - Gemergt wird erst nach Freigabe durch Alex.
2. **Dateibesitz.** Jedes Paket nennt die Dateien, die es ändern darf. Dateien eines
   anderen laufenden Pakets nicht anfassen. Wird das doch nötig: Arbeit anhalten, Alex fragen.
3. **Gemeinsame Dateien nur additiv.** `server.py`, `cli.py`, `units/__init__.py`,
   `web/src/types.ts`, `web/src/api.ts`, `README.md`: nur neue Einträge/Funktionen anhängen,
   nichts umbauen oder umformatieren. So bleiben Merge-Konflikte trivial.
4. **CHANGELOG.md** erst ganz am Ende des Pakets ergänzen (ein Eintrag unter `[Unreleased]`).
5. **Fertig heißt:** Tests für das Feature, `pytest` grün, bei Editor-Änderungen
   `npm test` und `npm run build` in `web/` grün, ein Beispiel in `examples/` wo sinnvoll,
   Massen-/Energiebilanz geprüft (`balances.py`).
6. **Nichts veröffentlichen:** kein `v*`-Tag, kein PyPI, Demo-Container bleibt aus.
7. **Fremde unfertige Änderungen** (uncommittete Dateien, Stashes) nicht committen,
   verwerfen oder überschreiben.

---

## Welle 5 · Vorrang für strukturelle Verbesserungen

Die zwölf Prozessstudien auf Fluxsheet 0.3.0 liefern die Ausgangsfälle. Im
Kondensations-Minimalfall der Apfelsaftstudie bewirkte ein vorgegebener Wärmestrom von
−1.498,766 kW bei einer Saccharosespur tatsächlich −2.160,069 kW, ohne eine Warnung zur
Vorgabe auszulösen. Weitere Befunde betreffen Stoffdaten-Ersatzannahmen, geschlossene
Kreisläufe und fehlende Modellzustände. Der sequenziell-modulare Solver bleibt die Basis.

### G9 – Druck-Enthalpie-Flash und Ergebnisprüfung · **ChatGPT** · Status: gemergt

1. Den Reinstoff- und Spurenfall aus der Apfelsaftstudie in eine kleine Regression
   übernehmen (`studies/2026-11-24/campaign-10/05-verdampfung/kondensation-spurprobe.json`).
   Den P-H-Flash an der Phasengrenze korrigieren und jeden zurückgegebenen P-H-Zustand
   gegen die geforderte Enthalpie prüfen. Kann der Zielzustand nicht bestimmt werden,
   den Fehlschlag mit Sollwert, Istwert und Residuum melden.
2. Bei fest vorgegebener Wärmeleistung am `Heater` und `HeatExchanger` die Vorgabe
   unabhängig vom berechneten Wärmestrom prüfen. Die Energiebilanz zusätzlich mit den
   tatsächlich berechneten Strömen bewerten. Rechenerfolg, Spezifikationserfüllung und
   Bilanzstatus in CLI, API und Editor getrennt ausgeben. Bestehende Ergebnisfelder
   und erfolgreiche Aufrufe kompatibel halten.

**Abnahme.** Der Minimalfall trifft in beiden Zweigen den geforderten Wärmestrom
innerhalb der numerischen Toleranz und liefert einen plausiblen Zweiphasenzustand.
Ein absichtlich injizierter Spezifikationsfehler bleibt auch bei geschlossener Bilanz
sichtbar. Das Fünf-Effekt-Modell wird erneut geprüft, und eventuelle verbleibende
Abweichungen werden für jeden betroffenen Wärmetauscher ausgewiesen.

- Dateien: `src/fluxsheet/thermo_package.py`, `src/fluxsheet/units/heat.py`,
  `src/fluxsheet/balances.py`, `src/fluxsheet/solver.py`, Ergebnisdarstellung in
  CLI, API und Editor, zugehörige Regressionstests.

### G10 – Stoffdaten je Fließbild und sichtbare Ersatzannahmen · **ChatGPT** · Status: offen

Binäre Literaturparameter und UNIFAC-Gruppenzuordnungen im Fließbild deklarieren können.
Bei einer Ergänzung der Gruppen gelten die vorhandenen Reinstoffkorrelationen weiter.
Quelle, Modellkonvention und Gültigkeitsbereich mitführen. Vor dem Lauf anzeigen, wo
Parameter fehlen, geschätzt oder idealisiert werden. Einen optionalen strengen Lauf für
Modelle anbieten, in denen solche Ersatzannahmen ausdrücklich ausgeschlossen werden.

**Abnahme.** Eine Fluxsheet-Datei enthält reproduzierbar dieselben Parameter in CLI
und Editor. Alte Dateien laden weiterhin. Für Ethanol/Glycerol und Wasser/Glycerol
den Vorrang lokaler Parameter, die Diagnose ohne diese Parameter und die Vorhersage
an separaten ternären Messpunkten prüfen.

- Dateien: `src/fluxsheet/flowsheet.py`, `src/fluxsheet/thermo_package.py`,
  `src/fluxsheet/fitting/store.py`, Stoffmethoden-Einstellungen im Editor und Tests.

### G11 – Stationäre geschlossene Kreisläufe · **ChatGPT** · Status: offen

Eine unabhängige Umlaufspezifikation auf einer internen Leitung zulassen. Den
numerischen Startwert für einen aufgetrennten Rückführstrom separat angeben können.
Ein geschlossener Kreis darf aus dieser Spezifikation starten. Bei einem unterbestimmten
Kreis erklären, welche Vorgabe fehlt.

**Abnahme.** Der R744-Referenzkreis rechnet ohne Zulauf bei spezifizierten 0,10 kg/s
Umlauf und reproduziert unter gleichen Annahmen Ströme und Wärmeleistungen des
bilanzierten Schnittmodells. Ein Modell ohne Umlaufspezifikation erhält eine gezielte
Diagnose. Bestehende Rückführungsmodelle mit Zulauf laufen weiterhin.

- Dateien: `src/fluxsheet/flowsheet.py`, `src/fluxsheet/solver.py`,
  Spezifikationseingabe im Editor und Rückführungs-Regressionstests.

### G12 – Schnittstellen für zusätzliche Phasen und Apparatezustände · **ChatGPT** · Status: offen

Für Feststoffphasen, reaktive Elektrolyte, zeitabhängige Adsorberbeladungen und
Polymer-Kettenmomente festlegen, welche Größen ein Stoffstrom transportiert und welche
dem Apparat gehören. Bilanzgrößen, Stoffmethoden, Serialisierung und Plug-in-Anbindung
dafür dokumentieren. Den Austausch am kleinsten Feststoff-Referenzfall von der
Thermodynamik über einen Strom bis zur Bilanz praktisch erproben.

**Abnahme.** Der Referenzfall unterscheidet Flüssigkeit und Feststoff, gibt beide
Phasen aus und bilanziert Masse und Enthalpie einschließlich Phasenwechsel. Bestehende
Dampf-Flüssig- und Flüssig-Flüssig-Modelle bestehen ihre Regressionstests. Für jeden
weiteren Modellzustand liegt eine dokumentierte Schnittstelle mit einer messbaren
Abnahmegröße vor.

- Dateien: `docs/architecture.md`, `src/fluxsheet/stream.py`,
  Thermodynamik- und Bilanzschnittstellen sowie Tests für den Referenzfall.

---

## Welle 1

### C1 – PI-Regler auf Pumpendrehzahl stabilisieren · **Claude** · Status: gemergt
Bekannter Fehler: Ein PI-Regler, der die Pumpendrehzahl stellt, versagt bei einem
Sollwertsprung. Ursache finden (Kennlinien-Auswertung außerhalb des Bereichs, Anti-Windup,
Steifigkeit), beheben, Regressionstest mit Sollwertsprung.
- Dateien: `src/fluxsheet/dynamics/*`, `src/fluxsheet/curves.py`,
  `src/fluxsheet/units/pressure_machines.py`, `tests/test_dynamics.py`, `tests/test_curves_logic.py`

### C3 – Ergebnis-Overlays im Fließbild · **Claude** · Status: gemergt
Nach dem Lauf T, p, Mengenstrom (und wahlweise Zusammensetzung) direkt an Strömen und
Units im Editor anzeigen, umschaltbar, mit Einheiten aus `quantities.ts`. Warnungen
(Bilanz, Phasengrenzen) als Marker an der Unit.
- Dateien: `web/src/StreamEdge.tsx`, `web/src/UnitNode.tsx`, neue `web/src/overlays.ts`
  (+ Test), `web/src/App.tsx` (nur Schalter), CSS

### G1 – Reaktive Gleichgewichtsstufen in der Kolonne · **ChatGPT** · Status: gemergt
Neben den kinetischen reaktiven Stufen: Stufen, auf denen gewählte Reaktionen das chemische
Gleichgewicht erreichen (K(T) wie im `EquilibriumReactor`), in die MESH-Gleichungen in
`columns/stages.py` eingebunden. Validierung gegen einen Literaturfall (z. B.
Methylacetat- oder MTBE-Reaktivdestillation), Beispiel in `examples/`.
- Dateien: `src/fluxsheet/columns/stages.py`, `src/fluxsheet/units/columns.py`,
  `src/fluxsheet/reactions.py` (nur Lesen/Hilfsfunktionen), `tests/test_columns.py`,
  `tests/test_column_side_streams.py`

### G2 – Nicht-idealer Gibbs-Reaktor · **ChatGPT** (übernommen von Claude) · Status: gemergt
`GibbsReactor` minimiert heute nur für ideales Gas. Erweitern um Fugazitäten aus dem
gewählten Property-Paket (Zustandsgleichung) und eine flüssige Phase (Aktivitäts-
koeffizienten), wahlweise mehrphasig. Test gegen bekannte Gleichgewichte bei hohem Druck
(z. B. Ammoniak- oder Methanolsynthese).
- Dateien: `src/fluxsheet/units/reactors.py` (nur `GibbsReactor`),
  `src/fluxsheet/thermo_package.py` (nur neue Hilfsfunktionen), `tests/test_reactors.py`

### G1b – Konvergenz reaktiver Gleichgewichtskolonnen · **ChatGPT** · Status: gemergt
Befund aus dem Review von G1: Das Beispiel `methyl_acetate_equilibrium_column.json`
konvergiert bei Rücklaufverhältnis 1,9 und 2,5, aber ab 3 nicht mehr (Molenbruch-Basis
schon ab 2,5); ohne Reaktion konvergiert dieselbe Kolonne bis 4. Robuster machen, z. B.
Homotopie über die Umsätze (ln K bzw. Umsatz von 0 schrittweise hochfahren) oder
Fortsetzung aus der Lösung bei kleinerem Rücklauf, wenn der direkte Newton-Lauf scheitert.
Fertig, wenn das Beispiel bei Rücklauf 2, 3, 4 und 5 auf beiden Basen konvergiert
(Regressionstest), alle Bilanzen schließen und die bestehenden Kolonnentests grün bleiben.
- Dateien: `src/fluxsheet/columns/stages.py` (nur `_initial`, `_newton`, `solve_column`
  und neue Hilfsfunktionen), `tests/test_column_side_streams.py`

> Konfliktpotential Welle 1: keins in Dateien. G1 und G2 nacheinander abarbeiten
> (erst G1 mergen), da beide `reactions.py` lesen könnten. G1b und C2 laufen parallel:
> C2 liest `columns/stages.py` nur und ändert die Datei nicht, solange G1b offen ist.

---

## Welle 2

### C2 – Dynamische Kolonne · **Claude** · Status: gemergt
Bodenhold-ups (Francis-Wehr), Kondensator-/Sumpffüllstand, Rückfluss- und Sumpfregelung
im bestehenden DAE-Kern; Start aus der stationären MESH-Lösung. Beispiel: Füllstands-
und Temperaturregelung einer Benzol/Toluol-Kolonne mit Feedsprung.
- Dateien: `src/fluxsheet/dynamics/*`, `src/fluxsheet/columns/stages.py` (nur Lesen bzw.
  Hilfsfunktionen herausziehen erst nach G1b),
  `src/fluxsheet/units/columns.py` (nur `dyn_params`), `tests/test_dynamics.py`,
  `web/src/DynamicsEditor.tsx`, `web/src/DynamicsView.tsx`

### G3 – Kostenrechnung für den Batch-Reaktor · **ChatGPT** (übernommen von Claude) · Status: gemergt
Turton-Kosten für den `BatchReactor` (Rührkessel, Mantel, Rührer) aus Volumen/Material;
Betriebskosten pro Charge und pro Jahr (Chargen/Jahr, Rüstzeit), damit `tac` auch für
Batch-Fließbilder gilt.
- Dateien: `src/fluxsheet/costing.py`, `src/fluxsheet/units/batch.py` (nur Kostendaten),
  `tests/test_optimization.py` bzw. neue `tests/test_costing_batch.py`

### G4 – DWSIM-Import · **ChatGPT** (übernommen von Claude) · Status: gemergt
Import von DWSIM-Dateien (`.dwxmz`/`.dwxml`): Komponenten, Property-Methode, Feeds,
Units, Verbindungen, Layout auf Fluxsheet-JSON abbilden; nicht abbildbare Teile als
Warnliste ausgeben. CLI `fluxsheet import-dwsim`, optional `POST /api/import`.
- Dateien: neues Paket `src/fluxsheet/importers/`, `src/fluxsheet/cli.py` (additiv),
  `src/fluxsheet/server.py` (additiv), `tests/test_import_dwsim.py`, Testdateien in `tests/data/`

---

## Welle 3

### C4 – Plug-ins über Entry-Points · **Claude** · Status: gemergt · *zuletzt, allein*
Fremde Pakete können Units, Property-Methoden und Komponenten über
`[project.entry-points."fluxsheet.units"]` usw. registrieren; Editor zeigt sie in der
Palette. Fasst die Registry an, daher als letztes Paket ohne parallele Unit-Arbeit.
- Dateien: `src/fluxsheet/units/__init__.py`, `src/fluxsheet/thermo_package.py`
  (Registry), `pyproject.toml`, `web/src/Palette.tsx`, `docs/architecture.md`,
  Beispiel-Plug-in unter `examples/plugin/`

### G5 – Validierungskampagne für die neuen Features · **ChatGPT** (übernommen von Claude) · Status: gemergt
G1, G2, G3, C2 gegen Referenzdaten in die bestehende Kampagne aufnehmen
(`docs/validation-campaign.json`, `docs/model-validation.md`), Abweichungen dokumentieren.
- Dateien: `docs/validation-campaign*`, `docs/model-validation.md`,
  `tests/test_validation_campaign.py`, Referenzdaten unter `src/fluxsheet/data/`

---

## Welle 4 (ab 2026-11-23)

Mehrere Agenten parallel: höchstens zwei Sessions gleichzeitig plus die koordinierende
Session. C6 folgt auf C5 und G7 auf G6, weil sie dieselben Dateien anfassen.

### C5 – Druckgetriebene Dynamik · **Claude** (Session) · Status: gemergt
Durchflüsse in der dynamischen Simulation folgen aus Druckdifferenzen statt fest
vorgegeben zu sein: Ventile mit Kennlinie (Cv, Öffnung als Stellgröße), Gas-Hold-up mit
Druck in Behältern (Druckaufbau, Abblasen), Pumpen gegen den Netzdruck. Eine Option
`pressure_driven` hält den bisherigen flussgetriebenen Modus als Standard.
- Dateien: `src/fluxsheet/dynamics/` (Kern, neue Module), `src/fluxsheet/units/pressure.py`
  (nur neue Ventil-Parameter, additiv), `tests/test_dynamics_pressure.py`, ein Beispiel in
  `examples/`, Abschnitt in `docs/architecture.md`

### C6 – Plug-ins überall · **Claude** (Session) · Status: gemergt
Plug-in-Units in der Dynamik (eigene dynamische Modelle über Entry-Point) und in der
Kostenrechnung; Henry-Komponenten mit Plug-in-Aktivitätsmethoden.
- Dateien: `src/fluxsheet/plugins.py`, `src/fluxsheet/dynamics/units.py` (Registry),
  `src/fluxsheet/costing.py` (Registry), `src/fluxsheet/thermo_package.py` (Henry-Pfad),
  `tests/test_plugins.py`, `examples/plugin/`

### G6 – Gekühlter Rohrreaktor · **ChatGPT** (übernommen von Claude, Session) · Status: gemergt
PFR mit Kühl- oder Heizmedium entlang der Länge im Gleich- oder Gegenstrom
(Randwertproblem beim Gegenstrom), Temperaturprofil beider Seiten, Prozessströme als
Mantelmedium (zweiter Strom am Reaktor).
- Dateien: `src/fluxsheet/units/reactors.py` (nur PFR), `src/fluxsheet/utilities.py`,
  `tests/test_pfr_cooled.py`, ein Beispiel in `examples/`

### G7 – Wärmedurchgang aus Filmkoeffizienten · **Claude** (übernommen von ChatGPT, Session) · Status: gemergt
U aus Geometrie und Korrelationen (Rohr innen: Gnielinski/Dittus-Boelter, Mantel, Rührkessel,
Kondensation/Sieden vereinfacht) plus Foulingwiderstände, statt U als feste Zahl.
- Dateien: `src/fluxsheet/utilities.py`, `src/fluxsheet/units/heat_transfer.py`,
  `src/fluxsheet/units/heat.py` (nur U-Herkunft), U-Eingabe des PFR in
  `src/fluxsheet/units/reactors.py`, neu `src/fluxsheet/film_coefficients.py` und
  `src/fluxsheet/film_rating.py`, `tests/test_film_coefficients.py`,
  `examples/double_pipe_from_geometry.json`; additiv: `thermo_package.py`
  (`transport_properties`), drei Größen in `quantities.py`, `report.py`,
  `web/src/quantities.ts`

### G8 – Restpunkte G1b und G4 · **ChatGPT** (übernommen von Claude) · Status: gemergt
Reaktive Gleichgewichtskolonnen mit Boilup-Spezifikation konvergieren lassen und die
Rücklauf-Fortsetzung beschleunigen (früher Abbruch, weniger Schritte); Import-Test mit einer
echten DWSIM-Datei, sofern eine mit passender Lizenz verfügbar ist.
- Dateien: `src/fluxsheet/columns/stages.py`, `src/fluxsheet/units/columns.py`,
  `src/fluxsheet/importers/dwsim.py`, `tests/test_column_side_streams.py`,
  `tests/test_import_dwsim.py`, `tests/data/dwsim/`

---

## Offen bei Alex (kein Agent)

- Drei ChatGPT-Stashes (`chatgpt-wip-before-…` zweimal, `Preserve local roadmap…`) sind
  nicht vollständig gemergt: behalten, anwenden oder löschen?
- Release 0.3.0 ist vorbereitet (Version, CHANGELOG); Tag, GitHub-Release, PyPI und
  Online-Demo erst auf ausdrücklichen Wunsch. 0.2.0 wurde nie veröffentlicht.

## Übersicht

| Welle | Claude | ChatGPT |
|---|---|---|
| 1 | C1 PI-Pumpenregler, C3 Ergebnis-Overlays | G1 Reaktive Gleichgewichtsstufen, G1b Konvergenz reaktiver Kolonnen, G2 Nicht-idealer Gibbs |
| 2 | C2 Dynamische Kolonne | G3 Batch-Kosten, G4 DWSIM-Import |
| 3 | C4 Plug-ins | G5 Validierungskampagne |
| 4 | C5 Druckgetriebene Dynamik, C6 Plug-ins überall | G6 Gekühlter Rohrreaktor, G7 Filmkoeffizienten, G8 Restpunkte G1b/G4 |
| 5 · als Nächstes | | G9 Flash und Ergebnisprüfung → G10 Stoffdaten → G11 geschlossene Kreisläufe → G12 Modellzustände |
