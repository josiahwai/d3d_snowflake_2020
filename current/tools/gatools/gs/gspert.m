%   USAGE: [response,eqx,misc] = gspert(eq,tok_data_struct,options,idoplots)
% 
%   PURPOSE: Calculate a plasma response with constraints on current
%            versus normalized area and thermal energy versus normalized volume
% 
%   INPUTS:
%     eq = structure containing equilibrium information
%     tok_data_struct = toksys structure containing name, geometry, greens.
%     options = optional structure with the following optional fields:
%        options.iconstraints (default = 1) selects a set of parameters that will
%        independently affect the plasma. The choices are (Is=conductor currents):
% 	  1. betap, li, Ip, Is.
% 	  2. Wth, li, Ip, Is.
%          3. Wth(V/Vtot), I(A/Atot), Is
%        If options is not structure it will be interpreted as options.iconstraints
%        The behavior of profiles can be chosen with the following:
%        options.dIndbetap, dIndw = recommended scaling of I vs A/Atot with betap or Wth
%        options.dIndli = recommended scaling of I versus A/Atot with li
%        options.dIndip = recommended scaling of I versus A/Atot with Ip
%        options.dWndbetap, dWndw = recommended scaling of W vs V/Vtot with betap or Wth
%        options.dWndli = recommended scaling of W versus V/Vtot with li
%        options.dWndip = recommended scaling of W versus V/Vtot with Ip
%        Returned response differs if for instance dWndbetap renders dli~=0, dbetap~=1, dip~=0
%        options.q = q values for which to calculate contours and their responses
%        options.gapspec(ngap,1:6), gap specification = (r,z,gr,gz,rs,zs),
%          i.e. a point (r,z), a vector (gr,gz), and an optional separatrix point (rs,zs)
%          A gap is *defined* as length of vector (gr,gz) extending from (r,z) to separatrix
%          If gapspec=(0,0,0,0,rs,zs), (gr,gz) is set to grad psi and (r,z) to a limiter point.
%        options.iso(niso,1:2), isoflux point specification = (r,z)
%        options.iverbose = flag to select text messages:
%          2=verbosinator, 1=talk, 0=shut up (default), -1=don't even warn
%        options.itdinclude = flag to include TD (=3D) in the group S of conductors
%          the default is 0 which sets ntd=0
%        options.idxcc = indices of CC elements in eigenmode analysis (default=1:ncc)
%        options.idxvv = indices of VV elements in eigenmode analysis (default=1:nvv)
%        options.idxtd = indices of TD (3D) elements in eigenmode analysis (default=1:ntd)
%        options.iwait = error handling: 0 (default) = warn, 1 = wait with message
%        options.iframe4movie = flag to make jpegs of idoplots, 0 = don't (default)
%        options.gspertfigno = figure number for plots, (111 is the default)
%        options.what2plot = choose what to show in the contour plots:
%          1. Change of psizr (default)
%          2. Change of normalized psi, dpsi-psibar*dpsibry-(1-psibar)*dpsimag
%          3. Change of flux generated by plasma currents, dpsizr_pla
%          4. Change of current density on grid, dcphi/dr/dz
%        options.nkn = number of knots for spline functions. Overrides the automatic number.
%     idoplots = (optional) array of flags to plot responses to individual quantities.
%        if idoplots is scalar, the value will apply to all plots: default is 0=don't plot
%        When iconstraints=1 or 2 the indices of idoplots correspond to:
%        1=w or betap, 2=li, 3=Ip, 4=time (unstable mode),
%        4+(1:nss)=conductors, 4+nss+(1:nss)=eigen modes
%        When iconstraints=3 the indices of idoplots correspond to:
%        [1:nr] = W profile, nr+[1:nr] =I profile, 2*nr+4=time (unstable mode),
%        2*nr+4+(1:nss)=conductors, 4+nss+(1:nss)=eigen modes
%        0 means no plot, <0 means wait, >0 is time in seconds to pause.
% 
%   OUTPUTS:
%     response = structure with responses.
%       Responses of the following quantities are provided:
%         current in grid cells (cphi),
%         current centroid (r,z),
%         boundary points (rb,zb) with touch or x-point as first element in the arrays,
%         all four possible strike points
%         axis and boundary flux, (psimag,psibry)
%         the functions pprime, ffprim,
%         the q-profile (qpsi),
%         contours of q = options.q (rq(nq,:), zq(nq,:)),
%         gap responses, how distance from (r,z) along (gr,gz) to separatrix changes
%         iso flux, br, bz responses
%         the geometric quantities L, A, V, I, W (explained below)
%       Here, the change of cphi, dcphi = sum(djphi*dr*dz)+dIedge
%       The responses are due to changes in these quantities:
%         coil and vessel currents (is),
%         total thermal energy (w) or betap or points on the profile W(V/Vtot),
%         li and total plasma current (ip) or points on the profile of I(A/Atot)
%       The efit definition of betap, li is used [Lao et al., Nucl. Fusion volume 25 (1985)]
%       In the case of qpsi and q-contours the response to Bphi (at eq.rzero) is added.
%       Response to time is also provided by:
%         Portone stability margin, growth rate and current vector (dIs) for unstable mode
%       points = indices to points on profiles of W(v), I(a) that are used as contraints
%     eqx = structure with extra equilibrium quantities derived from the equilibrium data:
%       (rcur, zcur), (ra, za), (r0, z0), (rx, zx), (rstrike, zstrike),
%         loci of current centroid, mag. axis, boundary defining point, x-points, strike points
%       drsep, distance between x-point flux surfaces measured in the outboard midplane
%       ilimited, a flag indicating if the plasma is limited, otherwise it is diverted
%       Rcnt(1:na,1:nr), Zcnt(1:na,1:nr), flux contours 
%       L(1:nr), surface integral of 1/R enclosed by flux contours
%       A(1:nr), area enclosed by flux contours
%       V(1:nr), volume enclosed by flux contours
%       I(1:nr), current enclosed by flux contours
%       W(1:nr), thermal energy enclosed by flux contours
%       betap = 4/3*mu0*W(end)/V(end)/bp2flx (should agree well with efit value)
%       li, (should agree well with efit value)
%       q = options.q (a direct copy of the input options.q)
%       psiq, the flux values [Wb] for the q values
%       qprime, dq/dpsi [1/(Wb/rad)]
%       Rq, Zq coordinates for contours of q
%       gaps(ngap,1:5) = [gap, r, z, rs, zs] where gap is distance from (r,z) on limiter to
%         (rs,zs) on the separatrix along vector [gr,gz] = [rs,zs]-[r,z];
%     misc = structure with miscellaneous information about the calculation:
%       execution_began, the date the code was called
%       texec, times at which parts were finished
%       mexec, message about what parts were finished
%       lots of others variables as well
% 
%   RESTRICTIONS:
% 
%   METHOD:
%     A perturbed Grad-Shafranov equation is solved subject to contraints
%     on the profiles of I(a), W(v),
%     where a, v are normalized area and volume within contours of poloidal flux
%     and I, W are plasma current and thermal energy within the contours.
%     The flux functions pres, and fpol^2/2 are allowed to change by 
%     spline functions (with d(pres)=0 and d(fpol^2/2) = R^2*Bphi*dBphi at edge)
%     A matrix of equations for the perturbed equilibrium is created and solved.
%     The equations are:
%       flux change on grid and current, thermal energy changes within contours
%     The variables being solved for are: 
%       flux change on grid and coeffs for spline functions
% 
%     Additional explanation:
%     The change of an axisymmetric equilibrium can be calculated exactly from the
%     Grad-Shafranov equation if changes in applied field, pprime, ffprim are known.
%     This code uses changes of thermal energy versus normalized volume, W(v) and
%     toroidal current versus normalized area, I(a) instead of pprime and ffprim.
%     An empirical scaling between the functions W(v), I(a) and a set of selectable
%     plasma parameters such as {betap, li, Ip} is utilized to completely define the
%     new equilibrium (with {dbetap, dli, dip, dis} as inputs in the example).
%     The returned plasma responses are such that they leave the other selectable
%     parameters unchanged, for instance, dcphidli does not change betap or Ip.
%     The set of selectable parameters is chosen with options.iconstraints where the
%     default 1 corresponds to {betap, li, Ip}
%