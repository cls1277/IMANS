function [v]=NDS(fit1,fit2) %������������֮���֧�����ϵ a֧��b ����1 b֧��a����2 �໥��֧�䷵��0
    v=0;
    dom_less=0;
    dom_equal=0;
    dom_more=0;
    for k=1:2
        if (fit1(1,k)>fit2(1,k))  
            dom_more = dom_more + 1;
        elseif (fit1(1,k)==fit2(1,k))
            dom_equal = dom_equal + 1;
        else
            dom_less = dom_less + 1;
        end
    end
    if dom_less == 0 && dom_equal ~= 2 % ˵��a��b֧��
        v=2;
    end
    if dom_more == 0 && dom_equal ~= 2 % ˵��a֧��b
        v=1;
    end
end