function  Xref = CorruptX(Xref, qqa,qqb,qqc,qq1,qq2,qq3)
    [a b] = size(Xref);
    switch nargin
        case 4
            for i = 1:a
            Xref(i,1) = Xref(i,1) + normrnd(0,sqrt(qqa));
            Xref(i,2) = Xref(i,2) + normrnd(0,sqrt(qqb)); 
            Xref(i,3) = Xref(i,3) + normrnd(0,sqrt(qqc));
            end

        case 7
            for i = 1:a
            Xref(i,1) = Xref(i,1) + normrnd(0,sqrt(qqa));
            Xref(i,2) = Xref(i,2) + normrnd(0,sqrt(qqb)); 
            Xref(i,3) = Xref(i,3) + normrnd(0,sqrt(qqc));
            Xref(i,4) = Xref(i,4) + normrnd(0,sqrt(qq1));
            Xref(i,5) = Xref(i,5) + normrnd(0,sqrt(qq2)); 
            Xref(i,6) = Xref(i,6) + normrnd(0,sqrt(qq3));
            end
        otherwise
            error('You need 4 or 7 input values!');
    end
end
