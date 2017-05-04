function[AICResC,BICResC,AICResKnot,BICResKnot]=Mix_VarCo_LLRR_PenBIC(Y,X,U,Z,V,Component,Knot,Int,lambda,penalty)
%Use BIC to determine the tuning parameter
%Component是选择的范围是1*C的一个向量
%首先获取自由度
[S,~,~]=size(V);
[P,~,~]=size(X);
[Q,~,~]=size(Z);
[~,SK]=size(Knot);
[~,SL]=size(lambda);
%其次获取样本量
N=sum(sum(Y~=inf));
AICResPenLoglikeli=-inf;
BICResPenLoglikeli=-inf;
for sk=1:SK
    temKnot=Knot(:,sk);
    KL=sum(temKnot)+4*P;
    for sl=1:SL
        temlambda=lambda(1,sl);
        [~,~,~,~,~,Loglikeli,temC]=Mixture_VaryingCo_LongReg_RanE_Sp(Y,X,U,Z,V,Component,temKnot,Int,temlambda,penalty);
        Degree=temC*(KL+Q+(S+1)*S/2+1+1)-1;
        temAICPenLoglikeli=Loglikeli-0.5*Degree*2;
        temBICPenLoglikeli=Loglikeli-0.5*Degree*log(N);
        if temAICPenLoglikeli>AICResPenLoglikeli
            AICReslambda=temlambda
            AICResC=temC
            AICResKnot=temKnot
            AICResPenLoglikeli=temAICPenLoglikeli
        end
        if temBICPenLoglikeli>BICResPenLoglikeli
            BICReslambda=temlambda
            BICResC=temC
            BICResKnot=temKnot
            BICResPenLoglikeli=temBICPenLoglikeli
        end
    end
end