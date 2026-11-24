import type { CheckStatus, SolveResult } from "./types";

const CHECK_LABELS: Record<CheckStatus, string> = {
  ok: "Passed",
  warning: "Out of tolerance",
  not_checkable: "Unable to check fully",
  not_applicable: "Not applicable",
  not_run: "Not run",
};

export function SolveStatus({ result }: { result: SolveResult }) {
  const { specification_status: specifications, balance_status: balances } = result;
  if (!specifications || !balances) return null;

  const incomplete = [specifications, balances].some(
    (status) => status === "warning" || status === "not_checkable",
  );
  return (
    <div className={`banner run-status${incomplete ? " warn" : ""}`} role="status">
      <span><b>Calculation</b> {result.converged ? "Converged" : "Not completed"}</span>
      <span><b>Duty specifications</b> {CHECK_LABELS[specifications]}</span>
      <span><b>Balances</b> {CHECK_LABELS[balances]}</span>
    </div>
  );
}
