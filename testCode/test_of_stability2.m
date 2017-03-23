

F_full = [G,zeros(5,110),-At;A,-eye(110,110),zeros(110,110);zeros(110,5),diag(lambda),diag(y)];

cond_F_full = cond(F_full)
[L,U] = lu(F_full);
cond_U = cond(U)

cond_F = cond(F)
L = chol(F);
cond_L = cond(L')