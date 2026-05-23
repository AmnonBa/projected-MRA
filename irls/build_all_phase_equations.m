function eqs = build_all_phase_equations(M3hat, rHat, p, clipCosines)
    q = (p - 1)/2;

    Blist = [];
    betaList = [];
    ampList = [];

    % Non-wrapping equations: a+b=c <= q
    for a = 1:q
        for b = a:q
            c = a + b;
            if c <= q
                amp = 2*rHat(a)*rHat(b)*rHat(c);
                if amp <= 1e-12
                    continue;
                end

                cosval = M3hat(a,b,c) / amp;
                if clipCosines
                    cosval = max(min(real(cosval),1),-1);
                elseif abs(cosval) > 1 || ~isfinite(cosval)
                    continue;
                end

                row = zeros(1,q);
                row(a) = row(a) + 1;
                row(b) = row(b) + 1;
                row(c) = row(c) - 1;

                Blist = [Blist; row]; 
                betaList = [betaList; acos(cosval)]; 
                ampList = [ampList; amp];
            end
        end
    end

    % Wrapping equations: a+b+c=p
    for a = 1:q
        for b = a:q
            for c = b:q
                if a + b + c == p
                    amp = 2*rHat(a)*rHat(b)*rHat(c);
                    if amp <= 1e-12
                        continue;
                    end

                    cosval = M3hat(a,b,c) / amp;
                    if clipCosines
                        cosval = max(min(real(cosval),1),-1);
                    elseif abs(cosval) > 1 || ~isfinite(cosval)
                        continue;
                    end

                    row = zeros(1,q);
                    row(a) = row(a) + 1;
                    row(b) = row(b) + 1;
                    row(c) = row(c) + 1;

                    Blist = [Blist; row]; 
                    betaList = [betaList; acos(cosval)]; 
                    ampList = [ampList; amp]; 
                end
            end
        end
    end

    eqs.B = sparse(Blist);
    eqs.beta = betaList;
    eqs.amp = ampList;
end
