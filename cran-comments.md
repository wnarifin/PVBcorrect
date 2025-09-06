## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission.

* Note:
"...no visible binding for '<<-' assignment to ..."
The use of <<- is required to use `weight_k` for weight update procedure.
The use of <<- is required to use `counter` to update bootstrap progress.
Otherwise, the functions will output wrong results. Other solutions were explored, but not good.
