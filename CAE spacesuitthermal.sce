
// ============================================================
// SPACESUIT LIQUID COOLING SYSTEM - THERMAL ANALYSIS
// ============================================================
//
// Purpose:
// Estimate heat transfer from the human body into a liquid
// cooled spacesuit garment using 16 meters of tubing.
//
// Thermal path:
// Body -> Tube/Body Gap -> Tube Wall -> Water
//
// Steady-state model:
// Q_body = Q_water
//
// Assumptions:
// - No radiation heat transfer
// - No external environmental heat loss
// - Uniform body temperature
// - Constant metabolic heat production
// - Constant water inlet temperature
// - Constant water properties
//
// Outputs:
// 1. Body temperature vs water flow rate
// 2. Coolant outlet temperature vs flow rate
// 3. Thermal resistance breakdown
// 4. Overall UA vs flow rate
// 5. CSV sensitivity results
//
// ============================================================


// ============================================================
// 1. INPUT PARAMETERS
// ============================================================

// Human body
Q_body = 100;                // Metabolic heat (W)

// Tubing geometry
L = 16;                      // Tube length (m)
D_outer = 0.008;             // Tube OD (m)
D_inner = 0.006;             // Tube ID (m)

wall_thickness = (D_outer-D_inner)/2;

// Tube material: polyurethane
k_tube = 0.20;               // W/m-K

// Coolant: water
rho_water = 997;             // kg/m3
cp_water = 4180;              // J/kg-K
mu_water = 0.00089;          // Pa-s
k_water = 0.60;               // W/m-K

Pr_water = cp_water*mu_water/k_water;

// Coolant inlet temperature
T_water_in = 15;             // C

// Nominal flow rate
flow_nominal_gpm = 0.5;

// Nominal gap conductance coefficient
h_gap_nominal = 50;          // W/m2-K

// Gap conductance sensitivity cases
gap_values = [10 25 50 100 200];


// ============================================================
// 2. GEOMETRY
// ============================================================

A_outer = %pi*D_outer*L;
A_inner = %pi*D_inner*L;

// Cylindrical tube wall conduction resistance
R_wall = log(D_outer/D_inner) / ...
         (2*%pi*k_tube*L);


// ============================================================
// 3. UNIT CONVERSION
// ============================================================

function flow_m3s = gpm_to_m3s(flow_gpm)

    flow_m3s = flow_gpm*3.785411784e-3/60;

endfunction


// ============================================================
// 4. INTERNAL CONVECTION CALCULATION
// ============================================================

function [Re,Nu,h_internal,velocity] = ...
         internal_convection(flow_gpm)

    // Convert flow to SI units
    volume_flow = gpm_to_m3s(flow_gpm);

    // Internal cross-sectional area
    area_flow = %pi*D_inner^2/4;

    // Average water velocity
    velocity = volume_flow/area_flow;

    // Reynolds number
    Re = rho_water*velocity*D_inner/mu_water;

    // Nusselt number
    if Re <= 2300 then

        // Laminar flow
        Nu = 3.66;

    elseif Re >= 3000 then

        // Turbulent flow
        Nu = 0.023*Re^0.8*Pr_water^0.4;

    else

        // Transition region
        Nu_laminar = 3.66;

        Nu_turbulent = ...
            0.023*3000^0.8*Pr_water^0.4;

        fraction = (Re-2300)/(3000-2300);

        Nu = Nu_laminar + ...
             fraction*(Nu_turbulent-Nu_laminar);

    end

    // Internal convection coefficient
    h_internal = Nu*k_water/D_inner;

endfunction


// ============================================================
// 5. THERMAL RESISTANCE NETWORK
// ============================================================

function [R_gap,R_wall_out,R_conv,R_total,UA,...
          Re,Nu,h_internal,velocity] = ...
          calculate_UA(flow_gpm,h_gap)

    // Internal convection
    [Re,Nu,h_internal,velocity] = ...
        internal_convection(flow_gpm);

    // Body-to-tube gap resistance
    R_gap = 1/(h_gap*A_outer);

    // Tube wall resistance
    R_wall_out = R_wall;

    // Water-side convection resistance
    R_conv = 1/(h_internal*A_inner);

    // Total thermal resistance
    R_total = R_gap+R_wall_out+R_conv;

    // Overall thermal conductance
    UA = 1/R_total;

endfunction


// ============================================================
// 6. SYSTEM HEAT TRANSFER CALCULATION
// ============================================================

function [T_body,T_out,Q_actual,UA,...
          R_gap,R_wall_out,R_conv,R_total,...
          Re,Nu,h_internal,velocity,mass_flow,...
          C_dot,effectiveness] = ...
          calculate_system(flow_gpm,h_gap)

    // Water volumetric flow
    volume_flow = gpm_to_m3s(flow_gpm);

    // Water mass flow
    mass_flow = rho_water*volume_flow;

    // Water heat capacity rate
    C_dot = mass_flow*cp_water;

    // Thermal resistance network
    [R_gap,R_wall_out,R_conv,R_total,UA,...
     Re,Nu,h_internal,velocity] = ...
        calculate_UA(flow_gpm,h_gap);

    // Heat exchanger effectiveness
    effectiveness = 1-exp(-UA/C_dot);

    // Assumed heat removal
    Q_actual = Q_body;

    // Required steady-state body temperature
    T_body = T_water_in + ...
             Q_actual/(C_dot*effectiveness);

    // Coolant outlet temperature
    T_out = T_water_in+Q_actual/C_dot;

endfunction


// ============================================================
// 7. NOMINAL OPERATING POINT
// ============================================================

[T_body,T_out,Q_actual,UA,...
 R_gap,R_wall_out,R_conv,R_total,...
 Re,Nu,h_internal,velocity,mass_flow,...
 C_dot,effectiveness] = ...
    calculate_system(flow_nominal_gpm,h_gap_nominal);


mprintf("\n");
mprintf("==========================================\n");
mprintf(" SPACESUIT COOLING SYSTEM - NOMINAL CASE\n");
mprintf("==========================================\n");

mprintf("Heat production:       %.2f W\n",Q_body);
mprintf("Tube length:           %.2f m\n",L);
mprintf("Tube OD:               %.2f mm\n",D_outer*1000);
mprintf("Tube ID:               %.2f mm\n",D_inner*1000);
mprintf("Water flow:            %.2f GPM\n",flow_nominal_gpm);
mprintf("Water inlet:           %.2f C\n",T_water_in);
mprintf("Gap conductance:       %.2f W/m2-K\n",h_gap_nominal);

mprintf("------------------------------------------\n");
mprintf("THERMAL RESISTANCES\n");

mprintf("Gap resistance:        %.5f K/W\n",R_gap);
mprintf("Tube wall resistance:  %.5f K/W\n",R_wall_out);
mprintf("Water convection:      %.5f K/W\n",R_conv);
mprintf("Total resistance:      %.5f K/W\n",R_total);

mprintf("------------------------------------------\n");
mprintf("FLOW CONDITIONS\n");

mprintf("Reynolds number:       %.0f\n",Re);
mprintf("Internal h:            %.2f W/m2-K\n",h_internal);
mprintf("Water velocity:        %.4f m/s\n",velocity);
mprintf("Water mass flow:       %.5f kg/s\n",mass_flow);
mprintf("Overall UA:            %.3f W/K\n",UA);

mprintf("------------------------------------------\n");
mprintf("THERMAL RESULTS\n");

mprintf("Heat removed:          %.2f W\n",Q_actual);
mprintf("Body temperature:      %.2f C\n",T_body);
mprintf("Coolant outlet:        %.2f C\n",T_out);

mprintf("==========================================\n");


// ============================================================
// 8. FLOW RATE SENSITIVITY SWEEP
// ============================================================

flow_rates = linspace(0.05,1.50,100);

n_flow = length(flow_rates);
n_gap = length(gap_values);

// Preallocate result matrices
body_temps = zeros(n_gap,n_flow);
outlet_temps = zeros(n_gap,n_flow);
UA_values = zeros(n_gap,n_flow);

for j = 1:n_gap

    h_gap = gap_values(j);

    for i = 1:n_flow

        flow = flow_rates(i);

        [Tb,To,Q,UA_temp,...
         Rg,Rw,Rc,Rt,...
         Re_temp,Nu_temp,hi,vel,mf,C,effect] = ...
            calculate_system(flow,h_gap);

        body_temps(j,i) = Tb;
        outlet_temps(j,i) = To;
        UA_values(j,i) = UA_temp;

    end

end


// ============================================================
// 9. GRAPH 1 - BODY TEMPERATURE VS FLOW RATE
// ============================================================

scf(1);
clf();

plot(flow_rates',body_temps');

xgrid();

xlabel("Water Flow Rate (GPM)");
ylabel("Required Body Temperature (deg C)");

title("Spacesuit Cooling: Body Temperature vs Water Flow");

legend( ...
    "10 W/m2-K", ...
    "25 W/m2-K", ...
    "50 W/m2-K", ...
    "100 W/m2-K", ...
    "200 W/m2-K", ...
    "in_upper_left");

a = gca();
a.font_size = 3;
a.thickness = 2;


// ============================================================
// 10. GRAPH 2 - COOLANT OUTLET TEMPERATURE
// ============================================================

scf(2);
clf();

plot(flow_rates',outlet_temps');

xgrid();

xlabel("Water Flow Rate (GPM)");
ylabel("Coolant Outlet Temperature (deg C)");

title("Coolant Outlet Temperature vs Flow Rate");

legend( ...
    "10 W/m2-K", ...
    "25 W/m2-K", ...
    "50 W/m2-K", ...
    "100 W/m2-K", ...
    "200 W/m2-K", ...
    "in_upper_left");

a = gca();
a.font_size = 3;
a.thickness = 2;


// ============================================================
// 11. GRAPH 3 - THERMAL RESISTANCE BREAKDOWN
// ============================================================

scf(3);
clf();

resistance_values = [R_gap R_wall_out R_conv];

bar(resistance_values);

xgrid();

xtitle( ...
    "Thermal Resistance Breakdown at 0.5 GPM", ...
    "Thermal Resistance Component", ...
    "Thermal Resistance (K/W)");

a = gca();
a.font_size = 3;


// ============================================================
// 12. GRAPH 4 - OVERALL UA VS FLOW RATE
// ============================================================

scf(4);
clf();

plot(flow_rates',UA_values');

xgrid();

xlabel("Water Flow Rate (GPM)");
ylabel("Overall UA (W/K)");

title("Overall Thermal Conductance vs Water Flow");

legend( ...
    "10 W/m2-K", ...
    "25 W/m2-K", ...
    "50 W/m2-K", ...
    "100 W/m2-K", ...
    "200 W/m2-K", ...
    "in_upper_left");

a = gca();
a.font_size = 3;
a.thickness = 2;


// ============================================================
// 13. EXPORT SENSITIVITY RESULTS TO CSV
// ============================================================

filename = "spacesuit_thermal_sensitivity.csv";

fid = mopen(filename,"wt");

mfprintf(fid, ...
    "Flow_GPM,Gap_Conductance_W_m2K," + ...
    "Body_Temperature_C,Coolant_Outlet_C,UA_W_K\n");

for j = 1:n_gap

    for i = 1:n_flow

        mfprintf(fid,"%.6f,%.3f,%.4f,%.4f,%.5f\n", ...
            flow_rates(i), ...
            gap_values(j), ...
            body_temps(j,i), ...
            outlet_temps(j,i), ...
            UA_values(j,i));

    end

end

mclose(fid);


// ============================================================
// 14. NOMINAL FLOW GAP CONDUCTANCE SENSITIVITY
// ============================================================

mprintf("\n");
mprintf("SENSITIVITY AT NOMINAL FLOW OF 0.5 GPM\n");
mprintf("------------------------------------------\n");

mprintf("Gap h       Body Temp       Outlet Temp\n");

for j = 1:n_gap

    [Tb,To,Q,UA_temp,...
     Rg,Rw,Rc,Rt,...
     Re_temp,Nu_temp,hi,vel,mf,C,effect] = ...
        calculate_system(flow_nominal_gpm,gap_values(j));

    mprintf("%6.1f      %8.2f C      %8.2f C\n", ...
        gap_values(j),Tb,To);

end

mprintf("\nAnalysis complete.\n");
mprintf("Four plots generated.\n");
mprintf("CSV file generated: %s\n",filename);

