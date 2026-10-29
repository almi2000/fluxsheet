// Symbols for the reactors on a 40x40 grid: a vessel with a mark for the model inside.
import { S } from "./style";

const vessel = <rect {...S} x="8" y="4" width="24" height="32" rx="7" />;
const thin = { ...S, strokeWidth: 1.5 };

export const shapes: Record<string, JSX.Element> = {
  StoichiometricReactor: (
    <>
      {vessel}
      <path {...thin} d="M13 20h14M23 16l4 4-4 4" />
    </>
  ),
  YieldReactor: (
    <>
      {vessel}
      <path {...thin} d="M14 29v-6M20 29v-12M26 29v-8" />
    </>
  ),
  EquilibriumReactor: (
    <>
      {vessel}
      <path {...thin} d="M13 17h14l-4-3M27 23H13l4 3" />
    </>
  ),
  GibbsReactor: (
    <>
      {vessel}
      <path {...thin} d="M13 12q7 22 14 0" />
      <circle cx="20" cy="22.5" r="1.8" fill="currentColor" />
    </>
  ),
  CSTR: (
    <>
      <rect {...S} x="7" y="9" width="26" height="27" rx="5" />
      <path {...thin} d="M20 3v22M13 25l7-3 7 3" />
    </>
  ),
  PFR: (
    <>
      <rect {...S} x="3" y="13" width="34" height="14" rx="7" />
      <path {...thin} d="M11 27l5-14M19 27l5-14M27 27l4-11" />
    </>
  ),
};
