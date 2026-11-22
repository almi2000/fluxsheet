# Fluxsheet – Roadmap und Arbeitsteilung (Claude / ChatGPT)

Stand: 2026-11-20 · Basis: `main` @ `85976bc`

Diese Datei verteilt die offenen Features auf zwei Agenten, die **gleichzeitig** arbeiten:
**Claude** (Claude Code) und **ChatGPT**. Die Pakete sind so geschnitten, dass in jeder
Welle die beiden Agenten möglichst **verschiedene Dateien** anfassen. Eine Welle ist fertig,
wenn beide Pakete in `main` gemergt sind; erst dann beginnt die nächste.

Status-Werte: `offen` · `in Arbeit` · `Review` · `gemergt`. Wer ein Paket beginnt oder
abschließt, ändert nur die Statuszeile **seines** Pakets.

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

### G2 – Nicht-idealer Gibbs-Reaktor · **ChatGPT** · Status: offen
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

### G3 – Kostenrechnung für den Batch-Reaktor · **ChatGPT** · Status: offen
Turton-Kosten für den `BatchReactor` (Rührkessel, Mantel, Rührer) aus Volumen/Material;
Betriebskosten pro Charge und pro Jahr (Chargen/Jahr, Rüstzeit), damit `tac` auch für
Batch-Fließbilder gilt.
- Dateien: `src/fluxsheet/costing.py`, `src/fluxsheet/units/batch.py` (nur Kostendaten),
  `tests/test_optimization.py` bzw. neue `tests/test_costing_batch.py`

### G4 – DWSIM-Import · **ChatGPT** · Status: Review
Import von DWSIM-Dateien (`.dwxmz`/`.dwxml`): Komponenten, Property-Methode, Feeds,
Units, Verbindungen, Layout auf Fluxsheet-JSON abbilden; nicht abbildbare Teile als
Warnliste ausgeben. CLI `fluxsheet import-dwsim`, optional `POST /api/import`.
- Dateien: neues Paket `src/fluxsheet/importers/`, `src/fluxsheet/cli.py` (additiv),
  `src/fluxsheet/server.py` (additiv), `tests/test_import_dwsim.py`, Testdateien in `tests/data/`

---

## Welle 3

### C4 – Plug-ins über Entry-Points · **Claude** · Status: offen · *zuletzt, allein*
Fremde Pakete können Units, Property-Methoden und Komponenten über
`[project.entry-points."fluxsheet.units"]` usw. registrieren; Editor zeigt sie in der
Palette. Fasst die Registry an, daher als letztes Paket ohne parallele Unit-Arbeit.
- Dateien: `src/fluxsheet/units/__init__.py`, `src/fluxsheet/thermo_package.py`
  (Registry), `pyproject.toml`, `web/src/Palette.tsx`, `docs/architecture.md`,
  Beispiel-Plug-in unter `examples/plugin/`

### G5 – Validierungskampagne für die neuen Features · **ChatGPT** · Status: offen
G1, G2, G3, C2 gegen Referenzdaten in die bestehende Kampagne aufnehmen
(`docs/validation-campaign.json`, `docs/model-validation.md`), Abweichungen dokumentieren.
- Dateien: `docs/validation-campaign*`, `docs/model-validation.md`,
  `tests/test_validation_campaign.py`, Referenzdaten unter `src/fluxsheet/data/`

---

## Offen bei Alex (kein Agent)

- Zwei ChatGPT-Stashes in `main` (`chatgpt-wip-before-…`, je 15 Dateien) sind nicht
  vollständig gemergt: behalten, anwenden oder löschen?
- 8 alte Worktrees unter `.claude/worktrees/` und gemergte `claude/*`-Branches aufräumen.
- Release 0.2.0 (Tag, GitHub-Release, PyPI) und Online-Demo: bewusst angehalten.

## Übersicht

| Welle | Claude | ChatGPT |
|---|---|---|
| 1 | C1 PI-Pumpenregler, C3 Ergebnis-Overlays | G1 Reaktive Gleichgewichtsstufen, G1b Konvergenz reaktiver Kolonnen, G2 Nicht-idealer Gibbs |
| 2 | C2 Dynamische Kolonne | G3 Batch-Kosten, G4 DWSIM-Import |
| 3 | C4 Plug-ins | G5 Validierungskampagne |
