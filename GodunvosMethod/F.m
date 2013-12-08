function [output] = F(ul,ur)
s = (ul+ur) / 2;
if ul >= ur
    if s > 0
        output = ul^2;
    else
        output = ur^2;
    end
else
    if ul > 0
        output = ul^2;
    elseif ur < 0
        output = ur^2;
    else
        output = 0;
    end
end