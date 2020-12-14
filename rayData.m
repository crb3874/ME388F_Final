classdef rayData
   properties
       
       %% Same ray used for both +/- along that direction
       
       %% -> inputs
       xs       %Initial x coordinate
       ys       %Initial y coordinate
       ox       % x component of direction
       oy       % y component of direction
       xint     % coordinates of intersections (x)
       yint     % coordinates of intersections (y)
       psi      % [psi1, psi2] end values of flux in plus and minus direc.
       v_trL    % Track lengths
       v_trI    % Need i,j indices of cells
       v_trJ    % that tracks are going through
       v_trA    % Track attenuation factors
       v_trTau  % Track tau factors
       
       %%%%%%%%% END PROPERTIES BLOCK %%%%%%%%%
       
   end
   
   methods
       
       function ins = rayData( xs , ys , ox , oy , ozV , trI0 , trJ0 , mD ) %mD meshData object
           
           % Initialize some data for handling tracks
           ins.xs           = xs;
           ins.ys           = ys;
           ins.ox           = ox;
           ins.oy           = oy;
           ins.xint         = [xs];
           ins.yint         = [ys];
           ins.v_trI        = [trI0];
           ins.v_trJ        = [trJ0];
           
           % Will iterate on these..
           idx              = trI0;
           jdx              = trJ0;
           xpos             = xs;
           ypos             = ys;
           
           t = 0;
           
           if ox > 0
                xbound = mD.x(end);
           elseif ox < 0
                xbound = mD.x(1);
           end
                
           while (sign(ox) * (xbound-xpos)) > 0 && (ypos < mD.y(end))
               
               t=t+1;
               ins.v_trI(t) = idx;
               ins.v_trJ(t) = jdx;
               
               % Get dx for current position.
               if ins.ox > 0
                   
                   % Going right ->
                   dx = abs(xpos - mD.x(idx+1));
                   
               elseif ins.ox < 0
                   
                   % <- Going left
                   dx = abs(xpos - mD.x(idx));
                   
               end
               
               % Get dy for current position
               dy = abs( dx * (ins.oy/ins.ox) );
               
               % Will need to compare this with distance to upper bound.
               
               dist_y = abs(ypos - mD.y(jdx+1));
               
               %% Three possible cases -- ray exits through side, through
               %% top, or through the corner exactly.

               if abs(dy - dist_y) < 1e-12
                   % Exits through corner (woah!) 
               
                   xpos = xpos + sign(ins.ox) * dx;
                   ypos = mD.y(jdx+1);
                   
                   idx = idx + sign(ins.ox);
                   jdx = jdx + 1;
                   
                   ins.v_trL = [ins.v_trL, abs(dx/ins.ox)];
                   
               elseif dy > dist_y
                   % Exits through top
                   xpos = xpos + dist_y * ins.ox/ins.oy;
                   ypos = mD.y(jdx+1);
                   
                   idx = idx;
                   jdx = jdx + 1;
                   
                   ins.v_trL = [ins.v_trL, abs(dist_y/ins.oy)];
                   
               elseif dy < dist_y
                   % Exits through side
                   xpos = xpos + sign(ins.ox) * dx;
                   ypos = ypos + dy;
                   
                   idx = idx + sign(ins.ox);
                   jdx = jdx;
                   
                   ins.v_trL = [ins.v_trL, abs(dx/ins.ox)];                
               end
               
               ins.xint = [ins.xint,xpos];
               ins.yint = [ins.yint,ypos];
               
           end
           
           %% Precompute some tau factors and attenuation factors.
           %% Each ray is seeded with all possible polar angle component,
           %% so we can properly project the track lengths above into
           %% the x-y plane for each oz. 
               
           ins.v_trTau  = zeros(length(ozV),t);    
           ins.v_trA    = zeros(length(ozV),t);
           
           for k = 1:t
               
               ix = ins.v_trI(k);
               jx = ins.v_trJ(k);

               ins.v_trTau(:,k) = abs( mD.st(ix,jx) * ins.v_trL(k) ./ ozV );
               ins.v_trA(:,k)   = exp( -ins.v_trTau(:,k) );
               
           end
           
       end
       
       %%%%%%%%% END METHODS BLOCK %%%%%%%%%
       
   end
end