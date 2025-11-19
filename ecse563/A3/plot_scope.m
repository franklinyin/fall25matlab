function plot_scope(out, scopeName, xlabel_str, ylabel_str, title_str, filename)
    fromSim = out.get(scopeName);
    figure;
    hold on; grid on;
    for i = 1:fromSim.numElements
        toPlot = fromSim.get(i);
        plot(toPlot.Values.Time, toPlot.Values.Data, 'DisplayName', toPlot.Name);
    end
    xlabel(xlabel_str); ylabel(ylabel_str);
    title(title_str); saveas(gcf, filename);
    legend('show');
end

