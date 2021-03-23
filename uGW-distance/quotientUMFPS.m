function [uXt, muXt] = quotientUMFPS(uX,muX,t)

uXt_big = uX;
uXt_big(uX<=t) = 0;

I = [1];
dI = uXt_big(1,:);

for k=1:length(uX)
    mx = max(dI);
    if mx == 0
        break
    end
    J = find(dI == mx);
    I = [I, J(1)];
    dI = min([dI;uXt_big(J(1),:)],[],1);
end


uXt = uX(I,I);
muXt = zeros(length(I),1);
for i = 1:length(I)
    temp = uX(i,:);
    K = temp<=t;
    muXt(i) = sum(muX(K));
end

