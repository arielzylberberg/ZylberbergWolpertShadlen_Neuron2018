classdef bayes_net_draw < handle
    % class to draw bayes nets (ariel.zylberberg@gmail.com)
    
    properties
        n_nodes = 0;
        radius = 0.2;
        
        nodes
        connections = []
        plates = []
        
        text_offset = 1.5 % multiplies radius
        
    end
    methods
        function obj = bayes_net_draw()
        end
        
        function add_node(obj,x,y,id,visible)
            if nargin<5 || isempty(visible)
                visible = 0;
            end
            r = obj.radius;
            h = circle(x,y,r);
            obj.n_nodes = obj.n_nodes+1;
            obj.nodes(obj.n_nodes).x = x;
            obj.nodes(obj.n_nodes).y = y;
            obj.nodes(obj.n_nodes).h = h;
            obj.nodes(obj.n_nodes).id = id;
            obj.nodes(obj.n_nodes).visible = visible; % default
            obj.nodes(obj.n_nodes).hyper = false;
            
        end
        
        
        function add_hyperparam(obj,x,y,id)
            r = obj.radius/2;
            h = hyper_box(x,y,r);
            obj.n_nodes = obj.n_nodes+1;
            obj.nodes(obj.n_nodes).x = x;
            obj.nodes(obj.n_nodes).y = y;
            obj.nodes(obj.n_nodes).h = h;
            obj.nodes(obj.n_nodes).id = id;
            obj.nodes(obj.n_nodes).visible = false; % default
            obj.nodes(obj.n_nodes).hyper = true;
            
        end
        
        
        function add_connection(obj,id_from,id_to)
            ids = {obj.nodes.id};
            i = find(ismember(ids,id_from));
            j = find(ismember(ids,id_to));
            x1 = obj.nodes(i).x;
            y1 = obj.nodes(i).y;
            x2 = obj.nodes(j).x;
            y2 = obj.nodes(j).y;
            hold on
            
            
            %plot([x1,x2],[y1,y2],'k');
            dx = x2-x1;
            dy = y2-y1;
            norm = sqrt(dx^2+dy^2);
            dx = dx/norm;
            dy = dy/norm;
            %x1 = x1+dx*0.2;
            %y1 = y1+dy*0.2;
            x2 = x2 - dx * 0.25;
            y2 = y2 - dy * 0.25;
            %h = plot([x1,x2],[y1,y2],'k');
            
            %h = quiver(x1,y1,x2-x1,y2-y1);
            h = arrow([x1,y1],[x2,y2],'Length',10,'TipAngle',20,'BaseAngle',50,...
                'facecolor',0.5*[1,1,1],'edgecolor',0.5*[1,1,1],'LineWidth',2);
            obj.connections(end+1) = h;
            
            
            
        end
        
        function format(obj)
            h = {obj.nodes.h};
            c1 = [148,187,228]/255;
            c2 = [0 125 197]/255;
            for i=1:length(h)
                if obj.nodes(i).hyper==1
                    set(h{i},'facecolor',c2,'edgecolor',c2);
                else
                    if obj.nodes(i).visible==1
                        set(h{i},'facecolor',c1,'edgecolor',c2);
                    else
                        set(h{i},'facecolor','w','edgecolor',c2);
                    end
                end
                uistack(h{i},'top');
                set(h{i},'LineWidth',2);
            end
            % send text upwards
            set([obj.nodes.string],'interpreter','latex','color','k','HorizontalAlignment','center','FontSize',15);
            uistack([obj.nodes.string],'top');
        end
        
        function add_text(obj,id,string,offset)
            if nargin<4 || isempty(offset)
                offset = obj.text_offset;
            end
                
            ids = {obj.nodes.id};
            i = find(ismember(ids,id));
            x = obj.nodes(i).x - obj.radius*offset;
            y = obj.nodes(i).y;
            obj.nodes(i).string = text(x,y,string);
            
        end
        
        function add_plate(obj,rect,offset,str)
            r = rect;
            r(1) = r(1)-offset;
            r(2) = r(2)-offset;
            r(3) = r(3)+2*offset;
            r(4) = r(4)+2*offset;
            h = rectangle('Position',r,'Curvature',[0,0]);
            obj.plates(end+1) = h;
            if ~isempty(str)
                text(r(1),r(2)+r(4),str,'HorizontalAlignment','left',...
                    'VerticalAlignment','top','FontSize',15,'interpreter','tex');
            end
        end
        
    end
end


function h = circle(x,y,r)
d = r*2;
px = x-r;
py = y-r;
h = rectangle('Position',[px py d d],'Curvature',[1,1]);
daspect([1,1,1])
end

function h = hyper_box(x,y,r)
d = r*2;
px = x-r;
py = y-r;
h = rectangle('Position',[px py d d],'Curvature',[.3,.3]);
end


