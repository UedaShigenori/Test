package DieselEngineLibrary
  package Components
    package Radiator
    
    
    model Cell1
    //table data
      parameter String tableName = "Tab1" "Table name on file or in function usertab (see docu)" annotation(
        Dialog(group = "Table data definition"));
      parameter String fileName = "C:\\Work\\2020\\DieselEngineLibrary\\csv\\data_K.txt" "File where matrix is stored" annotation(
        Dialog(group = "Table data definition", loadSelector(filter = "Text files (*.txt);;csv-files (*.csv)", caption = "Open file in which table is present")));
    
      parameter Modelica.Thermal.FluidHeatFlow.Media.Medium waterMedium = StreamConnectors.Media.LLC() "Cooling water" annotation(
        choicesAllMatching = true);
      parameter Modelica.Thermal.FluidHeatFlow.Media.Medium windMedium = StreamConnectors.Media.Air_75degC() "Cooling wind" annotation(
        choicesAllMatching = true);
      parameter DieselEngineLibrary.SIunits.MMLength Lrws "基準ラジエーター幅   : Lrws  (サンプルデータ)";
      parameter DieselEngineLibrary.SIunits.MMLength Lrw "ラジエータ幅";
    //output
      output Modelica.SIunits.ThermodynamicTemperature Twind_out;
      output Modelica.SIunits.ThermodynamicTemperature Twater_out;
    //instance
      StreamConnectors.Components.HeatedPipeOneDirection pipe_water(medium = waterMedium)  annotation(
        Placement(visible = true, transformation(origin = {0, 34}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
      StreamConnectors.Sensors.VFlow vFlow1(medium = waterMedium)  annotation(
        Placement(visible = true, transformation(origin = {0, 74}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
      StreamConnectors.Sensors.VFlow vFlow(medium = windMedium)  annotation(
        Placement(visible = true, transformation(origin = {-156, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.Convection convection annotation(
        Placement(visible = true, transformation(origin = {-32, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      StreamConnectors.Components.HeatedPipeOneDirection pipe_wind(medium = windMedium)  annotation(
        Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Interfaces.FluidPort wind_in annotation(
        Placement(visible = true, transformation(origin = {-200, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Interfaces.FluidPort wind_out annotation(
        Placement(visible = true, transformation(origin = {150, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Interfaces.FluidPort water_in annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Interfaces.FluidPort water_out annotation(
        Placement(visible = true, transformation(origin = {0, -102}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Gain gain(k = Lrw / Lrws * 1000) annotation(
        Placement(visible = true, transformation(origin = {-64, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Tables.CombiTable2D combiTable2D(fileName = fileName, tableName = tableName, tableOnFile = true) annotation(
        Placement(visible = true, transformation(origin = {-112, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(vFlow.port_b, pipe_wind.port_a) annotation(
        Line(points = {{-146, 0}, {-71, 0}}));
      connect(vFlow1.port_b, pipe_water.port_a) annotation(
        Line(points = {{0, 64}, {0, 45}}));
      connect(vFlow1.port_a, water_in) annotation(
        Line(points = {{0, 84}, {0, 98}}));
      connect(pipe_water.port_b, water_out) annotation(
        Line(points = {{0, 23}, {0, -102}}));
      connect(vFlow.port_a, wind_in) annotation(
        Line(points = {{-166, 0}, {-200, 0}}));
      connect(combiTable2D.y, gain.u) annotation(
        Line(points = {{-101, 52}, {-76, 52}}, color = {0, 0, 127}));
      connect(pipe_wind.port_b, wind_out) annotation(
        Line(points = {{-49, 0}, {-6, 0}, {-6, 6}, {6, 6}, {6, 0}, {150, 0}}));
      connect(convection.solid, pipe_wind.heatPort) annotation(
        Line(points = {{-42, 34}, {-60, 34}, {-60, 6}}, color = {191, 0, 0}));
      connect(gain.y, convection.Gc) annotation(
        Line(points = {{-52, 52}, {-32, 52}, {-32, 44}}, color = {0, 0, 127}));
      connect(pipe_water.heatPort, convection.fluid) annotation(
        Line(points = {{-6, 34}, {-22, 34}}, color = {191, 0, 0}));
//  Twind_out=wind_in.h_outflow/windMedium.cp;
//  Twater_out=water_in.h_outflow/waterMedium.cp;
      Twind_out = pipe_wind.port_b.h_outflow / windMedium.cp;
      Twater_out=pipe_water.port_b.h_outflow/waterMedium.cp;
      connect(vFlow.m, combiTable2D.u1) annotation(
        Line(points = {{-156, 10}, {-156, 10}, {-156, 58}, {-124, 58}, {-124, 58}}, color = {0, 0, 127}));
      connect(vFlow1.m, combiTable2D.u2) annotation(
        Line(points = {{-8, 74}, {-18, 74}, {-18, 90}, {-148, 90}, {-148, 46}, {-124, 46}, {-124, 46}}, color = {0, 0, 127}));
    
    annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {150, 100}}, initialScale = 0.1), graphics = {Text(origin = {-57, 73}, extent = {{-23, 11}, {33, -9}}, textString = "幅を変更する場合のKの見直しと"), Text(origin = {-113, 73}, extent = {{-23, 11}, {23, -11}}, textString = "基準サンプルの熱通過率"), Text(origin = {-57, 63}, extent = {{-23, 11}, {19, -1}}, textString = "kW->Wへの単位変換"), Text(origin = {-121, 55}, extent = {{-23, 11}, {-7, 1}}, textString = "冷却風量"), Text(origin = {-121, 43}, extent = {{-25, 13}, {-7, 1}}, textString = "冷却水流量"), Text(origin = {-185, -21}, extent = {{-23, 11}, {-7, 1}}, textString = "冷却風"), Text(origin = {15, 105}, extent = {{-25, 13}, {-7, 1}}, textString = "冷却水"), Line(origin = {-199, -20}, points = {{-7, 0}, {7, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {19.7635, 104.15}, points = {{-7, 0}, {-7, -14}}, arrow = {Arrow.None, Arrow.Filled})}),
        Icon(coordinateSystem(extent = {{-110, -100}, {100, 100}})));
    end Cell1;
    

    
    model HeatExchanger
      //parameter
      parameter DieselEngineLibrary.SIunits.MMLength Dy(min = 0) = 25 "Cell ⊿y ・・高さピッチ";
      parameter DieselEngineLibrary.SIunits.MMLength Dz(min = 0) = 2 "Cell ⊿z ・・厚みピッチ";  
      parameter DieselEngineLibrary.SIunits.MMLength Lrw(min = 0) = 700 "ラジエータ幅" annotation(Dialog(group = "搭載ラジエータ"));
      parameter DieselEngineLibrary.SIunits.MMLength Lrh(min = 0) = 700 "ラジエータ高さ" annotation(Dialog(group = "搭載ラジエータ"));
      parameter DieselEngineLibrary.SIunits.MMLength Lrd(min = 0) = 46 "ラジエータ厚さ" annotation(Dialog(group = "搭載ラジエータ"));
      parameter DieselEngineLibrary.SIunits.MMLength Lrws(min = 0) = 525 "ラジエータ厚さ" annotation(Dialog(group = "基準ラジエータ"));
    //port
      DieselEngineLibrary.Interfaces.FluidPort wind_in annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Interfaces.FluidPort water_in annotation(
        Placement(visible = true, transformation(origin = {0, 98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Interfaces.FluidPort water_out annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -102}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      DieselEngineLibrary.Interfaces.FluidPort wind_out annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      StreamConnectors.Components.Cell1 cell[Nj, Nk](each Lrw=Lrw,each Lrws=Lrws) annotation(
        Placement(visible = true, transformation(origin = {1.83479, 0.582059}, extent = {{-10.7634, -9.41794}, {8.07252, 9.41794}}, rotation = 0)));
    
    
    
    protected
    //  parameter Integer Nj(min = 2) = integer(Lrh/Dy) "Cell 分割数　高さ(j)方向の数 (搭載ラジエータ用)";
        //    parameter Integer Nk(min = 2) = integer(Lrd/Dz) "Cell 分割数　厚み(k)方向の数 (搭載ラジエータ用)";
        //debug mode
      parameter Integer Nj(min = 2) = 28 "Cell 分割数　高さ(j)方向の数 (搭載ラジエータ用)";
      parameter Integer Nk(min = 2) = 23 "Cell 分割数　厚み(k)方向の数 (搭載ラジエータ用)";
    equation
//wind_in,wind_outポートと端側のCellの接続
        for n in 1:Nj loop
          connect(wind_in, cell[n, 1].wind_in);
          connect(wind_out, cell[n, Nk].wind_out);
        end for;
//water_in,water_outポートと上下側のCellの接続
        for m in 1:Nk loop
          connect(water_in, cell[1, m].water_in);
          connect(water_out, cell[Nj, m].water_out);
        end for;
//Cell同士の接続
        for p in 1:Nj - 1 loop
          for q in 1:Nk loop
            connect(cell[p, q].water_out, cell[p + 1, q].water_in);
          end for;
        end for;
    
      for r in 1:Nk - 1 loop
        for s in 1:Nj loop
          connect(cell[s, r].wind_out, cell[s, r + 1].wind_in);
        end for;
      end for;
      
      
    annotation(Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}), Polygon(origin = {-14, -17}, fillColor = {173, 173, 173}, fillPattern = FillPattern.Solid, points = {{-60, 79}, {-60, -57}, {60, -79}, {60, 57}, {-60, 79}}), Polygon(origin = {59, -17}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Solid, points = {{-13, 57}, {13, 79}, {13, -61}, {-13, -79}, {-13, -13}, {-13, 57}}), Polygon(origin = {-1, 62}, fillColor = {218, 218, 218}, fillPattern = FillPattern.Solid, points = {{-73, 0}, {-47, 22}, {73, 0}, {47, -22}, {-73, 0}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14, 31}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, 8.275}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.15, -14.475}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.3, -37.2}, points = {{-60, 11}, {60, -11}}), Line(origin = {-14.275, -61.5}, points = {{-60, 11}, {60, -11}}), Line(origin = {59, 31}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.025, 8.625}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.85, -14.475}, points = {{-13, -11}, {13, 11}}), Line(origin = {58.875, -36.85}, points = {{-13, -11}, {13, 11}}), Line(origin = {59.075, -61.5}, points = {{-13, -11}, {13, 11}}), Line(origin = {6.825, 69}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-0.85, 63.425}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {-8, 57.35}, points = {{-58.825, 11}, {61.175, -11}}), Line(origin = {54, -21.775}, points = {{0, 69}, {0, -69}}), Line(origin = {60.675, -17.375}, points = {{0, 69}, {0, -69}}), Line(origin = {67.7, -11.75}, points = {{0, 69}, {0, -69}}), Line(origin = {-72.5086, 83.4913}, points = {{-13, -12}, {-1, -2}}, arrow = {Arrow.None, Arrow.Filled}), Text(origin = {-80, 87}, extent = {{-6, 5}, {6, -5}}, textString = "k"), Line(origin = {-85.5909, 57.5}, points = {{0, 14}, {0, 0}}, arrow = {Arrow.None, Arrow.Filled}), Text(origin = {-90, 67}, extent = {{-6, 5}, {6, -5}}, textString = "j"), Line(origin = {-79.5357, 69.2262}, points = {{-6, 2}, {6, -2}}, arrow = {Arrow.None, Arrow.Filled}), Text(origin = {-72, 71}, extent = {{-6, 5}, {6, -5}}, textString = "i")}, coordinateSystem(initialScale = 0.1)),
        preferredView = "text");
    end HeatExchanger;
    
    
    
    end Radiator;

    package Commons
      model HeatedPipeOneDirection "A simple static pipe with heating/cooling"
    
        parameter Modelica.Thermal.FluidHeatFlow.Media.Medium medium =  Modelica.Thermal.FluidHeatFlow.Media.Medium() "Medium in the component" annotation(
          choicesAllMatching = true);
        parameter Real K = dp_nominal / m_flow_nominal ^ 2 "Pressure drop coefficient";
        parameter Real dp_nominal = 0.5 "Nominal pressure drop [bar]" annotation(
          Dialog(group = "Nominal values"));
        parameter Real m_flow_nominal = 1 "Nominal mass flow rate [kg/s]" annotation(
          Dialog(group = "Nominal values"));
        DieselEngineLibrary.Interfaces.FluidPort port_a annotation(
          Placement(transformation(extent = {{-120, -10}, {-100, 10}}), iconTransformation(extent = {{-120, -10}, {-100, 10}})));
        DieselEngineLibrary.Interfaces.FluidPort port_b annotation(
          Placement(transformation(extent = {{100, -10}, {120, 10}}), iconTransformation(extent = {{100, -10}, {120, 10}})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation(
          Placement(visible = true, transformation(origin = {-2, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
// Mass balance
        port_a.m_flow + port_b.m_flow = 0;
// Momentum balance
        port_a.p - port_b.p = K * port_a.m_flow * abs(port_a.m_flow);
// Enthalpy
        port_b.h_outflow = inStream(port_a.h_outflow) + heatPort.Q_flow / port_a.m_flow;
        port_a.h_outflow = port_b.h_outflow;
//calculate heatPort.T
        inStream(port_a.h_outflow) = heatPort.T * medium.cp;
        annotation(
          preferredView = "text",
          Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Rectangle(lineColor = {28, 108, 200}, fillColor = {28, 108, 200}, fillPattern = FillPattern.Solid, extent = {{-100, 40}, {100, -40}}), Rectangle(fillColor = {215, 215, 215}, fillPattern = FillPattern.Backward, extent = {{-100, 60}, {100, 40}}), Rectangle(fillColor = {215, 215, 215}, fillPattern = FillPattern.Backward, extent = {{-100, -40}, {100, -60}}), Line(points = {{-40, -80}, {40, -80}}, thickness = 5), Line(origin = {5.72368, 0.197369}, points = {{30, -74}, {40, -80}, {30, -86}}, thickness = 5)}),
          Diagram(coordinateSystem(preserveAspectRatio = false)),
          Documentation(info = "<html>
      <p>This is probably the simplest model of a heated/cooled pipe.</p>
      <h4>Governing equations</h4>
      <p>The equations &quot;on paper&quot; are:</p>
      <p style=\"margin-left: 30px;\">m<sub>in</sub> = m<sub>out</sub></p>
      <p style=\"margin-left: 30px;\">p<sub>in</sub> - p<sub>out</sub> = deltaP=K<sup>.</sup>m<sub>in</sub>|m<sub>in</sub>|</p>
      <p style=\"margin-left: 30px;\">m<sub>out</sub><sup>.</sup>h<sub>out</sub> = m<sub>in</sub><sup>.</sup>h<sub>in</sub> + Q</p>
      <p>Representing the conservation of <i>mass</i>, <i>momentum</i> and <i>energy</i>.</p>
      <h4>Implementation</h4>
      <p>To understand the implementation and the concept of stream connectors:</p>
      <ol>
      <li>comment all code in the <code>equation</code> section.</li>
      <li>&apos;check&apos; the model (in Dymola, hit the <b>F8</b> key). The code checker will say that the model is unbalanced with <b>six</b> unknowns and <b>two</b> equations.</li>
      <li>try to understand why <b>four</b> equations are needed when the model can be described with only <b>three</b> equations on paper.</li>
      </ol>
      <p><br>The fact is that in Modelica, for each stream connector you need to provide a value for the outgoing stream variable (<code>h_outflow</code> in this case) when or if the mass flow leaves the connector. </p>
      <p>In the pipe case it means that you must provide an expression for <code>port_b.h_outflow</code> (obviously) and <code>port_a.h_outflow</code> (for the case of flow reversal).</p>
      <p>On the Github Wiki page on <a href=\"https://github.com/justnielsen/ModelicaTutorials/wiki/Stream-connectors\">Stream connectors</a>, I&apos;ve provided one more example similar to this.</p>
      </html>"));
      end HeatedPipeOneDirection;
    end Commons;
  end Components;

  package Media "Medium properties"
    extends Modelica.Icons.MaterialPropertiesPackage;

    record Air_40degC "Medium: properties of air at 40 degC and 1 bar"
      extends Modelica.Thermal.FluidHeatFlow.Media.Medium(rho = 1.091, cp = 1010, cv = 720, lamda = 0.0264, nue = 16.3E-6);
      annotation(
        defaultComponentPrefixes = "parameter",
        Documentation(info = "<html>
  Medium: properties of air at 40&deg;C and 1 bar
  </html>"));
    end Air_40degC;

    record Air_75degC "Medium: properties of air at 75 degC and 1 bar"
      extends Modelica.Thermal.FluidHeatFlow.Media.Medium(rho = 0.983, cp = 1007, cv = 720, lamda = 0.0264, nue = 16.3E-6);
      annotation(
        defaultComponentPrefixes = "parameter",
        Documentation(info = "<html>
    Medium: properties of air at 75&deg;C and 1 bar
    </html>"));
    end Air_75degC;

    record LLC "Medium: properties of LLC"
      extends Modelica.Thermal.FluidHeatFlow.Media.Medium(rho = 1035, cp = 3320, cv = 720, lamda = 0.0264, nue = 16.3E-6);
      annotation(
        defaultComponentPrefixes = "parameter",
        Documentation(info = "<html>
    Medium: properties of LLC
    </html>"));
    end LLC;
  end Media;

  package Examples
  extends Modelica.Icons.ExamplesPackage;
    model RadTest1
  StreamConnectors.Sources.VolumeFlowBoundary_T volumeFlowBoundary_T(T = 313.15, V_flow = 0.98 / 28 * 28, medium = air_75degC) annotation(
        Placement(visible = true, transformation(origin = {-68, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      StreamConnectors.Sources.VolumeFlowBoundary_T volumeFlowBoundary_T1(T = 373.15, V_flow = 1.7037037037037E-03 / 23 * 23, medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-10, 66}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      StreamConnectors.Sources.PressureBoundary_T pressureBoundary_T(medium = air_75degC) annotation(
        Placement(visible = true, transformation(origin = {42, 2}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      StreamConnectors.Sources.PressureBoundary_T pressureBoundary_T1(medium = llc) annotation(
        Placement(visible = true, transformation(origin = {-8, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      parameter StreamConnectors.Media.LLC llc annotation(
        Placement(visible = true, transformation(origin = {71, 77}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      parameter StreamConnectors.Media.Air_75degC air_75degC annotation(
        Placement(visible = true, transformation(origin = {71, 57}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      StreamConnectors.Components.HeatExchanger heatExchanger annotation(
        Placement(visible = true, transformation(origin = {-8, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(volumeFlowBoundary_T1.port, heatExchanger.water_in) annotation(
        Line(points = {{-10, 56}, {-8, 56}, {-8, 12}, {-8, 12}}));
      connect(heatExchanger.water_out, pressureBoundary_T1.port) annotation(
        Line(points = {{-8, -8}, {-8, -55}}));
      connect(heatExchanger.wind_out, pressureBoundary_T.port) annotation(
        Line(points = {{2, 2}, {31, 2}}));
      connect(volumeFlowBoundary_T.port, heatExchanger.wind_in) annotation(
        Line(points = {{-56, 2}, {-18, 2}, {-18, 2}, {-18, 2}}));
      annotation(
        Diagram(graphics = {Text(origin = {-31, 21}, extent = {{-25, 15}, {17, -9}}, textString = "車両番号 2相当のパラメータ使用"), Text(origin = {-31, 21}, extent = {{-25, 15}, {17, -9}}, textString = "車両番号 2相当のパラメータ使用"), Text(origin = {83, 83}, extent = {{-25, 15}, {5, 1}}, textString = "流体物性"), Rectangle(origin = {73, 71}, extent = {{-27, 29}, {27, -29}}), Text(origin = {-53, -21}, extent = {{-23, 11}, {-7, 1}}, textString = "冷却風"), Text(origin = {7, 73}, extent = {{-25, 13}, {-7, 1}}, textString = "冷却水"), Line(origin = {-51.2764, -2.82918}, points = {{-4.72361, -1.17082}, {5.27639, -1.17082}, {1.27639, 0.82918}}), Line(origin = {-48, -5}, points = {{2, 1}, {-2, -1}}), Line(origin = {-51.2764, -2.82918}, points = {{-4.72361, -1.17082}, {5.27639, -1.17082}, {1.27639, 0.82918}}), Line(origin = {-48, -5}, points = {{2, 1}, {-2, -1}}), Line(origin = {-9.57881, 47.4138}, rotation = -90, points = {{-4.72361, -1.17082}, {5.27639, -1.17082}, {1.27639, 0.82918}}), Line(origin = {-9.57881, 47.4138}, rotation = -90, points = {{-4.72361, -1.17082}, {5.27639, -1.17082}, {1.27639, 0.82918}}), Line(origin = {-11.8669, 44.2493}, rotation = -90, points = {{2, 1}, {-2, -1}})}, coordinateSystem(initialScale = 0.1)));
    
    end RadTest1;
  end Examples;

  
  package Interfaces
  extends Modelica.Icons.InterfacesPackage;
  
  connector FluidPort "Simple stream connector."
    Real p "Potential/effort variable";
    flow Real m_flow "Flow variable";
    stream Real h_outflow "Specific enthalpy";
  
    annotation (
      preferredView="info",
      defaultComponentName="port",
      Icon(coordinateSystem(preserveAspectRatio=false), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={28,108,200},
            fillPattern=FillPattern.Solid)}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
  <p>The simples possible stream connector using </p>
  <ul>
  <li><span style=\"color: #1c6cc8;\">mass flow rate</span> as &quot;flow&quot; variable (in a connection point &quot;flow&quot; variables are automatically summed to zero like Kirchhoff&apos;s current law)</li>
  <li><span style=\"color: #1c6cc8;\">pressure</span> as &quot;potential&quot; (or &quot;effort&quot;) variable (in a connection point the potentials are equal)</li>
  <li><span style=\"color: #1c6cc8;\">specific enthalpy</span> as &quot;stream&quot; variable (property carried with the direction of the flow)</li>
  </ul>
  <h4>Specific enthalpy vs. temperature</h4>
  <p>Why use specific enthalpy instead of temperature as stream variable when temperature seems more intuitive?</p>
  <p><br>The answer is &quot;automatic mixing equation&quot;: When you connect more than two stream connectors, the Modelica tool can automatically calculate the mixing enthalpy from knowledge of the flow directions, since m1*h1 + m2*h2 + ... + mN*hN = 0. This also means that it is <b>not strictly necessary</b> to create specific mixer/splitter components.</p>
  <h4>Why &apos;h_outflow&apos;</h4>
  <p>Why is the enthalpy named <code>h_outflow</code> and not just <code>h</code>?</p>
  <p>The reason is, that the variable is only relevant when the flow &apos;goes out&apos; of the component to which the connector is attached. </p>
  <p>Yes, really. If you want to know the value of the ingoing enthalpy after a simulation, you shold look at the value of <code>h_outflow</code> of the adjacent component (from which the mass flow &apos;goes out&apos;). </p>
  <p>If a component needs to access the ingoing enthalpy it should use one of the functions <code>inStream()</code> or <code>actualStream()</code>.</p>
  </html>"));
  
  end FluidPort;
  
  
  end Interfaces;

  package Sensors
    extends Modelica.Icons.SensorsPackage;

    model VFlow "Reads Volume flow rate"
  parameter Modelica.Thermal.FluidHeatFlow.Media.Medium medium = Modelica.Thermal.FluidHeatFlow.Media.Medium() "Medium in the component" annotation(
        choicesAllMatching = true);
      Interfaces.FluidPort port_a annotation(
        Placement(transformation(extent = {{-120, -10}, {-100, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
      Interfaces.FluidPort port_b annotation(
        Placement(transformation(extent = {{100, -10}, {120, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
      Modelica.Blocks.Interfaces.RealOutput m annotation(
        Placement(transformation(extent = {{-216, -2}, {-196, 18}}), iconTransformation(extent = {{-10, -10}, {10, 10}}, rotation = 90, origin = {0, 90})));
    equation
// Balance equations
      port_a.m_flow + port_b.m_flow = 0;
      port_a.p = port_b.p;
      port_a.h_outflow = inStream(port_b.h_outflow);
      port_b.h_outflow = inStream(port_a.h_outflow);
      m = port_a.m_flow / medium.rho;
      annotation(
        Icon(coordinateSystem(preserveAspectRatio = false), graphics = {Line(points = {{0, 80}, {0, 20}}, color = {0, 0, 0}, pattern = LinePattern.Dot), Line(points = {{-100, 0}, {100, 0}}, color = {0, 0, 0}), Line(points = {{-60, 40}, {-40, 20}, {42, 20}, {60, 40}}, color = {0, 0, 0}), Line(points = {{-60, -40}, {-40, -20}, {42, -20}, {60, -40}}, color = {0, 0, 0})}),
        Diagram(coordinateSystem(preserveAspectRatio = false)));
    end VFlow;
  end Sensors;

  package SIunits
  type MMLength = Real(final quantity = "Length", final unit = "mm");

  end SIunits;
  annotation(
    uses(Modelica(version = "3.2.3")));
end DieselEngineLibrary;
