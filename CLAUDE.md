@AGENTS.md

## Claude Code in this repository

- Work in a worktree under `.claude/worktrees/` on a branch `claude/<package-id>`. ChatGPT
  works in the main checkout on `chatgpt/<package-id>`, so the main checkout is usually not
  on `main`.
- To merge into `main` without touching that checkout: merge in a temporary worktree on a
  detached `main`, run the checks there, then move `main` only if it has not changed:
  `git update-ref refs/heads/main <new> <old>`.
- At most two parallel sessions at a time.
- Run the full pytest suite in the background; it takes about 8 minutes.
