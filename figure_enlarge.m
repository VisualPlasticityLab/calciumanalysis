function h=figure_enlarge(h,mag)
    p=h.Position;
    set(h, 'Position', [p(1)-p(3)*(mag-1) p(2)-p(4)*(mag-1) p(3)*mag p(4)*mag])
end