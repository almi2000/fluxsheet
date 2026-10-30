import { BaseEdge, getSmoothStepPath, useStore, type EdgeProps } from "@xyflow/react";
import { memo } from "react";

const STUB = 20; // horizontal run out of an outlet / into an inlet
const CLEARANCE = 48; // distance between a recycle lane and the units above it
const R = 8; // corner radius

/** Orthogonal polyline through `pts` with rounded corners. */
function roundedPath(pts: [number, number][]): string {
  let d = `M ${pts[0][0]},${pts[0][1]}`;
  for (let i = 1; i < pts.length - 1; i++) {
    const [px, py] = pts[i - 1];
    const [x, y] = pts[i];
    const [nx, ny] = pts[i + 1];
    const r = Math.min(R, Math.hypot(x - px, y - py) / 2, Math.hypot(nx - x, ny - y) / 2);
    const ax = x - Math.sign(x - px) * r;
    const ay = y - Math.sign(y - py) * r;
    const bx = x + Math.sign(nx - x) * r;
    const by = y + Math.sign(ny - y) * r;
    d += ` L ${ax},${ay} Q ${x},${y} ${bx},${by}`;
  }
  const [lx, ly] = pts[pts.length - 1];
  return `${d} L ${lx},${ly}`;
}

/**
 * A stream between ports. Streams that run forward are smooth-step edges; streams that run
 * backward (recycles) leave to the right, drop below every unit they pass and come back in
 * from the left, so they never cross over the units on the main line.
 */
function StreamEdgeView(props: EdgeProps) {
  const { sourceX, sourceY, targetX, targetY, sourcePosition, targetPosition, markerEnd, style } = props;
  const labelProps = {
    label: props.label,
    labelStyle: props.labelStyle,
    labelShowBg: props.labelShowBg,
    labelBgStyle: props.labelBgStyle,
    labelBgPadding: props.labelBgPadding,
    labelBgBorderRadius: props.labelBgBorderRadius,
  };
  const backward = targetX < sourceX - STUB;
  const left = Math.min(sourceX, targetX) - STUB;
  const right = Math.max(sourceX, targetX) + STUB;
  // Lowest unit bottom between the two ends (a number, so the edge re-renders only when it moves).
  const floor = useStore((s) => {
    if (!backward) return 0;
    let max = Math.max(sourceY, targetY);
    for (const n of s.nodeLookup.values()) {
      const { x, y } = n.internals.positionAbsolute;
      const w = n.measured.width ?? 0;
      if (x + w >= left && x <= right) max = Math.max(max, y + (n.measured.height ?? 0));
    }
    return max;
  });

  if (!backward) {
    const [path, labelX, labelY] = getSmoothStepPath({ sourceX, sourceY, targetX, targetY, sourcePosition, targetPosition });
    return <BaseEdge path={path} labelX={labelX} labelY={labelY} {...labelProps} markerEnd={markerEnd} style={style} />;
  }
  const lane = floor + CLEARANCE;
  const path = roundedPath([
    [sourceX, sourceY],
    [sourceX + STUB, sourceY],
    [sourceX + STUB, lane],
    [targetX - STUB, lane],
    [targetX - STUB, targetY],
    [targetX, targetY],
  ]);
  return (
    <BaseEdge path={path} labelX={(sourceX + targetX) / 2} labelY={lane} {...labelProps} markerEnd={markerEnd} style={style} />
  );
}

export const StreamEdge = memo(StreamEdgeView);
