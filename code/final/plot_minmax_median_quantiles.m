function [medianplot_array, new_YTick, YTickLabel] = plot_minmax_median_quantiles(lineStyle,mins,maxs,medians,quant25s,quant75s,choose_logy,method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict)

      medianplot_array = zeros(1, length(method_array));  % for legend
      num_iter_array = 1:iter_save:maxiter;
      num_methods = length(method_array);
      
      hold on

      for i = 1:num_methods
          
          minmaxcolor_i = minmaxcolor_dict(method_array{i});

          if choose_logy
            h = fill([num_iter_array  fliplr(num_iter_array)], [log10(maxs(:,i)')  fliplr(log10(mins(:,i))')], minmaxcolor_i,'EdgeColor', 'none');
            set(h,'facealpha', .5)
            quantcolor_i = quantcolor_dict(method_array{i});
            h = fill([num_iter_array  fliplr(num_iter_array)], [log10(quant75s(:,i)')  fliplr(log10(quant25s(:,i))')], quantcolor_i,'EdgeColor', 'none');
            set(h,'facealpha', .5)
            medianplot_array(i) = plot( num_iter_array, log10(medians(:,i)),linecolor_dict(method_array{i}),'LineWidth',2,...
                       'DisplayName',displayname_dict(method_array{i}), 'LineStyle', lineStyle );
            ylabel('(log scale)')
          else 
            h = fill([num_iter_array  fliplr(num_iter_array)], [maxs(:,i)'  fliplr(mins(:,i)')], minmaxcolor_i,'EdgeColor', 'none');
            set(h,'facealpha', .5)
            quantcolor_i = quantcolor_dict(method_array{i});
            h = fill([num_iter_array  fliplr(num_iter_array)], [quant75s(:,i)'  fliplr(quant25s(:,i)')], quantcolor_i,'EdgeColor', 'none');
            set(h,'facealpha', .5)
            medianplot_array(i) = plot( num_iter_array,medians(:,i),linecolor_dict(method_array{i}),'LineWidth',2,...
                       'DisplayName',displayname_dict(method_array{i}), 'LineStyle', lineStyle );          
          end

      end
      
      % on y axis: replace t by 10^t for interesting t values    
      
      min_tick = floor(log10(min(min(mins))));
      max_tick = ceil(log10(max(max(maxs))));
      new_YTick = [0, min_tick, round(0.5*(min_tick+max_tick))];                       
      new_YTick = sort(unique(new_YTick));
      for ii = 1:length(new_YTick)
          if new_YTick(ii) == 0
              YTickLabel{ii} = '0';
          else
              YTickLabel{ii} = num2str(new_YTick(ii), '10^{%d}');
          end
      end
      
    ylim([min_tick max_tick])
    set(gca, 'YTick', new_YTick, 'YTickLabel', YTickLabel);
    drawnow
      
    axis square
      
    % remove legend
    leg = legend('figure()');
    set(leg,'visible','off')

    hold off
      
    end