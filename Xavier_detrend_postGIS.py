import openpyxl as xl
import arcpy
import os
from os import *
from arcpy import env
from openpyxl import Workbook
from openpyxl import load_workbook
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import csv

###### INPUTS ######
# excel file containing xyz data for station points
direct = r'Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO2\COMID17573013'
xyz_table = direct + '\\XY_elevation_table_100_smooth_3_spaced.xlsx'
centerline = direct + '\\las_files\\centerline\\smooth_centerline.shp'
station_lines = direct + '\\las_files\\centerline_sp4ft_sm100ft\\smooth_centerline_XS_4x250ft.shp'
DEM = direct + '\\las_files\\ls_nodt.tif'
process_footprint = direct + '\\las_footprint.shp'
detrend_workplace = direct + '\\LINEAR_DETREND_BP1960_3ft_spacing'
spatial_ref = arcpy.Describe(process_footprint).spatialReference
######

#Import xl file
wb = load_workbook(xyz_table)
ws = wb.active
print("Workbook " + str(xyz_table) + " loaded")

#Fill lists with necessary data
id = []
location = []
x = []
y = []
z = []
listoflist = [location, id, z, x, y]
listofcolumn = ["D", "A", "L", "I", "J"]

for i in range(0, len(listofcolumn)):
    for cell in ws[listofcolumn[i]]:
        listoflist[i].append(cell.value)
    del listoflist[i][0]

print(listoflist)

#define station point spacing, number of station points and index
point_spacing = int(location[1]- location[0])
number_of_points = int(location[-1]/point_spacing)

print("Point spacing: " + str(point_spacing))
print("Number of points: " + str(number_of_points))

location_np = np.array(location)
z_np = np.array(z)

def quadratic_fit(location_np, location, z_np, ws):
    print("Applying quadratic fit...")
    constants = np.polyfit(location_np, z_np, 2)

    fit_params = []
    for i in constants:
        fit_params.append(i)
    print(fit_params)

    #Use the found regresssion down channel
    fitted_z = []
    for i in location:
        z_fit_temp = (fit_params[0]*i**(2) + fit_params[1]*i + fit_params[2])
        fitted_z.append(z_fit_temp)

    print(fitted_z)

    residual = []
    if len(fitted_z) == len(z):
        print("Lists are same length, continue as planned...")
        for i in range(len(fitted_z)):
            residual.append(z[i] - fitted_z[i])

    else:
        print("Lists are wrong lengths/not matching, something is wrong...")

    mean_z = (sum(z) / len(z))
    squared_real = []
    squared_res = []

    #Calculate the total sum of squares, sum of residuals, and R^2
    for i in range(len(residual)):
        squared_real.append((z[i] - mean_z)**2)
        squared_res.append(residual[i]**2)

    print("List of squared variance from mean: " + str(squared_real))
    print("List of squared residuals: " + str(squared_res))

    R_squared = 1 - (sum(squared_res)/sum(squared_real))
    print("The coefficient of determination is: " + str(R_squared))

    # Calculate mean and SD for residual to calculate residual z scores as a list
    mean_res = sum(residual)/len(residual)
    sd_res = np.std(residual)
    res_z_score = [(z_res_value - mean_res) * 1.0 / sd_res for z_res_value in residual]
    print(residual)
    print("The mean residual is %.2f" %mean_res)

    # Add fitted z values to the xyz table
    cell_test = ws["D1"]
    print(cell_test.value)
    cell_test.value = "Quadratic fit z values"
    print(cell_test.value)

    if ws["D1"].value == "Quadratic fit z values":
        for i in range(2, len(fitted_z)):
            cell = ws.cell(row=i, column=4)
            cell.value = float(fitted_z[i])
    else:
        print("Something is wrong with writing to the excel sheet")


    wb.save(filename=xyz_table)

    print("Excel file ready for Arc processing!")

def linear_fit(location, z, ws, list_of_breakpoints=[]):
    # Applies a linear fit to piecewise sections of the longitudinal profile, each piece is stored in split_list
    # ADD 0 BEFORE ANY ADDED BREAKPOINTS OR THE FUNCTION WILL FAIL!!!!!!!!!!
    print("Applying linear fit")
    list_of_breakpoints.append(location[-1])
    print(list_of_breakpoints)
    fit_params = []
    split_location_list = []
    split_z_list = []
    # Split by breakpoints into a list of lists
    if len(list_of_breakpoints) > 0:
        slope_break_indices = [int(distance/point_spacing) for distance in list_of_breakpoints]
        for i in range(1, len(slope_break_indices)):
            temp_location_list = []
            temp_z_list = []
            temp_break_index = slope_break_indices[i]
            if i == 1:
                temp_break_index = slope_break_indices[i]
                split_location_list.append(location[:temp_break_index])
                split_z_list.append(z[:temp_break_index])
            #elif i < (len(slope_break_indices)):   FIX THIS SO WE CAN WORK WITH MORE THAN JUST ONE BREAK POINT
                #temp_last_break_index = slope_break_indices[i-1]
                # Splice locations into split_location_list
                #temp_location_list.append(location[temp_last_break_index:temp_break_index])
                #split_location_list.append(temp_location_list)
            else:
                print("In the else: " + str(temp_z_list))
                split_location_list.append(location[slope_break_indices[i-1]:])
                split_z_list.append(z[slope_break_indices[i-1]:])

        print("Breakpoints added...")
        print("Slope break index: " + str(slope_break_indices))
        print("Split location list: " + str(split_location_list))
        print("Split z list: " + str(split_z_list))

        # Get fit parameters for each section of the data
        if len(split_z_list) == len(split_location_list):
            for i in range(len(split_location_list)):
                temp_loc_array = np.array(split_location_list[i])
                temp_z_array = np.array(split_z_list[i])
                m, b = np.polyfit(temp_loc_array, temp_z_array, 1)
                fit_params.append([m, b])
            print("Fit param list: " + str(fit_params))
        else:
            print("Something went wrong, list lengths do not match...")

    else:
        m, b = np.polyfit(location, z, 1)
        fit_params.append([m, b])
        print("Fit params [m, b]: " + str(fit_params))

    # Make list of lists storing fitted z values
    list_of_lengths = []
    z_fit_list = []
    for i in range(len(split_z_list)):
        list_of_lengths.append(len(split_z_list[i]))

    # Add fitted z's into a list
    print(list_of_lengths)
    i = 0
    while i < len(list_of_lengths):
        print(i)
        print(list_of_lengths[i])
        for j in range(list_of_lengths[i]):
            z_fit_list.append(split_location_list[i][j]*fit_params[i][0]+fit_params[i][1])
        i += 1

    if len(z_fit_list) == len(z):
        print("List lengths are compatible, next we export to excel")
    else:
        print("Something went wrong, length of z =/ z_fit_list")
    print(z_fit_list)

    #Calculate residual and R^2
    residual = []
    for i in range(len(z_fit_list)):
        residual.append(z[i]-z_fit_list[i])
    print(residual)
    mean_z = (sum(z) / len(z))
    squared_real = []
    squared_res = []

    #Calculate the total sum of squares, sum of residuals, and R^2
    for i in range(len(residual)):
        squared_real.append((z[i] - mean_z)**2)
        squared_res.append(residual[i]**2)

    print("List of squared variance from mean: " + str(squared_real))
    print("List of squared residuals: " + str(squared_res))

    R_squared = 1 - (sum(squared_res)/sum(squared_real))
    print("The coefficient of determination is: " + str(R_squared))

    # Calculate mean and SD for residual to calculate residual z scores as a list
    mean_res = sum(residual)/len(residual)
    sd_res = np.std(residual)
    res_z_score = [(z_res_value - mean_res) * 1.0 / sd_res for z_res_value in residual]
    print(residual)
    print("The mean residual is %.2f" %mean_res)

    # Add fitted z values to the xyz table
    '''NOTE: ADJUST SHEET ROW FOR CELL_TEST AND COLUMN FOR THE CELL = WS.CELL COMMAND'''
    cell_test = ws["M1"]
    print(cell_test.value)
    cell_test.value = ("z_fit_%s" % (list_of_breakpoints[1]))
    print(cell_test.value)

    if ws["M1"].value == cell_test.value:
        print("Sheet activated...")
        for i in range(2, len(z_fit_list)):
            cell = ws.cell(row=i, column=13)
            cell.value = float(z_fit_list[i])
    else:
        print("Something is wrong with writing to the excel sheet")
    ws["A1"].value == "OID"


    wb.save(filename=xyz_table)

    print("Excel file ready for Arc processing!")

    return [fit_params, z_fit_list, residual, R_squared]

##### Detrend in arcgis #####
def detrend_that_raster(fit_z_xl_file, dem, footprint, spatial_ref, list_of_breakpoints=[]):
    # Turn fitted xl file to a csv
    wb = xl.load_workbook(fit_z_xl_file)
    ws = wb.active
    csv_file = 'Stationpoints_xyz_sp4ft.csv'
    arcpy.env.extent = footprint
    arcpy.env.workspace = direct
    arcpy.env.snapRaster = dem
    arcpy.overwriteoutput = True

    if not os.path.exists(process_footprint):
        os.makedirs(process_footprint)

    with open(csv_file, 'w', newline="") as f:
        col = csv.writer(f)
        for row in ws.rows:
            col.writerow([cell.value for cell in row])

    for i in list_of_breakpoints:
        column = ('z_fit_%s' % i)
        detrended_raster_file = detrend_workplace + "\\ras_detren.tif"
        points = arcpy.MakeXYEventLayer_management(csv_file, "POINT_X", "POINT_Y", out_layer=("fitted_station_points%s" % i), spatial_reference=spatial_ref, in_z_field=column)
        points = arcpy.SaveToLayerFile_management(points, ("fitted_station_points%s_sp4ft" % i).replace('.csv', '.lyr'))
        points = arcpy.CopyFeatures_management(points)
        print("Creating Thiessen polygons...")
    # Delete non-relevent fields tp reduce errors
        fields = [f.name for f in arcpy.ListFields(points)]
        dont_delete_fields = ['FID', 'Shape', 'POINT_X', 'POINT_Y', column]
        fields2delete = list(set(fields) - set(dont_delete_fields))
        points = arcpy.DeleteField_management(points, fields2delete)

        thiessen = arcpy.CreateThiessenPolygons_analysis(points, "thiespoly_sp4ft.shp", fields_to_copy='ALL')
        z_fit_raster = arcpy.PolygonToRaster_conversion(thiessen, column, ('theis_raster_%s_sp4ft.tif' % i), cell_assignment="CELL_CENTER")
        detrended_DEM = arcpy.Raster(dem) - arcpy.Raster(z_fit_raster)
        detrended_DEM.save(detrended_raster_file)
        print("DEM DETRENDED!")


###### Define plotting function ######
def diagnostic_quick_plot(location_np, z_np):
    x_plot = location_np
    y_plot = z_np
    plt.plot(x_plot, y_plot, 'r', label="Actual elevation profile")
    plt.xlabel("Thalweg distance downstream (ft)")
    plt.ylabel("Bed elevation (ft)")
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    return plt.show()

def make_quadratic_fit_plot(location_np, z_np, fit_params):
    # Make a plot of the quadratic fit
    # Plot longitudinal elevation profile and detrending fit line
    x_plot = location_np
    y1_plot = z_np
    y2_plot = (fit_params[0]*x_plot**(2) + fit_params[1]*x_plot + fit_params[2])
    plt.plot(x_plot, y1_plot, 'r', label="Actual elevation profile")
    plt.xlabel("Thalweg distance downstream (ft)")
    plt.ylabel("Bed elevation (ft)")
    plt.title("Quadratic detrended elevation profile")
    plt.plot(x_plot, y2_plot, 'b', label= "Elevation profile detrending quadratic: %s * x^2 + %.4f * x + %.4f" % (fit_params[0], fit_params[1], fit_params[2]))
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.legend(loc=1)
    return plt.show()

def make_linear_fit_plot(location_np, z_np, fit_params):
    x_plot = location_np
    y1_plot = z_np
    y2_plots = []
    for list in fit_params:
        y2_plots.append(list[0]*x_plot + list[1])
    plt.plot(x_plot, y1_plot, 'r', label="Actual elevation profile")
    plt.xlabel("Thalweg distance downstream (ft)")
    plt.ylabel("Bed elevation (ft)")
    plt.title("Linear piecewise detrended elevation profile")
    for i in range(len(y2_plots)):
        plt.plot(x_plot, y2_plots[i], 'b', label="Elevation profile linear detrending: %.4f * x + %.4f" % (y2_plots[i][0], y2_plots[i][1]))
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.legend(loc=1)
    return plt.show()

def make_residual_plot(location_np, residual, R_squared):
    '''Plot residuals across longitudinal profile, show R^2. Inputs are a numpy array of location, a list of residuals, and a float for R-squared'''
    x_plot = location_np
    y_plot = np.array(residual)
    y_zero = 0*x_plot
    plt.scatter(x_plot, y_plot, s=1, c='r')
    plt.plot(x_plot, y_zero, c='b')
    plt.xlim(0, max(location_np))
    plt.ylim(min(y_plot), max(y_plot))
    plt.xlabel("Distance down stream centerline (ft)")
    plt.ylabel("Residual")
    plt.title("Residuals: R^2 = %.4f" % R_squared)
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    return plt.show()


################## CALL FUNCTIONS AS NECESSARY ####################
#make_quadratic_fit_plot(location_np, z_np, fit_params)
#make_quadratic_residual_plot(location_np, residual, R_squared)
#quadratic_fit(location_np, location, z_np, ws)

diagnostic_quick_plot(location_np, z_np)
#fit_params = linear_fit(location, z, ws, list_of_breakpoints=[0, 3200])[0]
#residual = linear_fit(location, z, ws, list_of_breakpoints=[0, 3200])[2]
#R_squared = linear_fit(location, z, ws, list_of_breakpoints=[0, 3200])[3]
#make_residual_plot(location_np, residual, R_squared)
#make_linear_fit_plot(location_np, z_np, fit_params)
#detrend_that_raster(xyz_table, DEM, process_footprint, spatial_ref, list_of_breakpoints=[3200])

####### WHEN WE RETURN FIGURE OUT HOW TO TURN THIS INTO SOMETHING WE CAN DETREND THE DEM WITH ######