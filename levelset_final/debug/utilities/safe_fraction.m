function value = safe_fraction(num, den)
%SAFE_FRACTION Divide safely and return NaN when denominator is nonpositive.

    if den <= 0
        value = NaN;
    else
        value = num / den;
    end
end
