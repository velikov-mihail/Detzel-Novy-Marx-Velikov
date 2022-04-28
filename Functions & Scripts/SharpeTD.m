function Sharpe=SharpeTD(wF,rets,FacdWs,tcosts)
%Take weight in factors and make Sharpe ratio net of costs accounting
%for cost diversification.
%     wF=round(wF,5);
%     wF=wF./sum(wF);
    grossrets = rets*wF; %[Tx1]s

    wFnonMKT = wF(2:end);
    
    N=size(tcosts,2);
    portdW = FacdWs*kron(wFnonMKT,speye(N)); 
    PortTCsCD=nansum(abs(portdW).*tcosts,2); %Tx1; Equation (12)

    netrets = grossrets-PortTCsCD;

    Sharpe = mean(netrets)/std(netrets);
end 