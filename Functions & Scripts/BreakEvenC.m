function [c_scalar,success]=BreakEvenC(gross_factors1,tc1,gross_factors2,tc2,dates,startend,eps)
    if nargin < 7
        eps=0.00001;
    end

    s=startend(1);
    e=startend(2);
    
    options = optimoptions('fmincon','Display','off');

    fun=@(c_scalar_)(SqSharpeDiff(c_scalar_, gross_factors1, tc1, gross_factors2, tc2,[s e]));
    x0 = 1;
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [0];
    ub = [1];
    nonlcon = [];
    [c_scalar,SRDiff] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

    if SRDiff<eps
        success=1;
    else
        success=0;
        c_scalar=0;
    end
end