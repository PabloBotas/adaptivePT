#include "energy_range.hpp"

// Energies MeV
std::vector<float> energy = {
		0.000, 4.500, 5.000, 5.500, 6.000, 6.500, 7.000,
		7.500, 8.000, 8.500, 9.000, 9.500, 10.00, 12.50,
		15.00, 17.50, 20.00, 22.50, 25.00, 27.50, 30.00,
		35.00, 40.00, 45.00, 50.00, 55.00, 60.00, 65.00,
		70.00, 75.00, 80.00, 85.00, 90.00, 95.00, 100.0,
		105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0,
		140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0,
		175.0, 180.0, 185.0, 190.0, 195.0, 200.0, 205.0,
		210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0,
		245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0,
		280.0, 285.0, 290.0, 295.0, 300.0, 350.0, 400.0
};
// Janni82 range mm water
std::vector<float> rangeJanni = {
		0.0000, 0.3059, 0.3676, 0.4342, 0.5058, 0.5822, 0.6634,
		0.7493, 0.8399, 0.9351, 1.0348, 1.1391, 1.2479, 1.8575,
		2.5739, 3.3939, 4.3147, 5.3336, 6.4487, 7.6579, 8.9593,
		11.833, 15.057, 18.620, 22.515, 26.730, 31.259, 36.093,
		41.226, 46.650, 52.360, 58.348, 64.610, 71.138, 77.929,
		84.960, 92.263, 99.802, 107.57, 115.64, 123.84, 132.33,
		141.08, 150.03, 159.19, 168.49, 178.20, 187.92, 197.92,
		208.07, 218.44, 228.99, 239.74, 250.72, 261.85, 273.18,
		284.65, 296.34, 308.07, 320.13, 332.21, 344.54, 357.00,
		369.73, 382.57, 395.41, 408.58, 421.82, 435.30, 448.86,
		462.64, 476.44, 490.41, 504.42, 518.71, 668.12, 829.02
};

Energy_Range_Calculator_t::Energy_Range_Calculator_t(Direction_t dir_) : dir(dir_)
{
    table.append("range", rangeJanni);
    table.append("energy", energy);
}


float Energy_Range_Calculator_t::calculate(float x, Direction_t dir)
{
    switch(dir)
    {
        case Energy_Range_Calculator_t::FromEtoR:
            return table.getVal("range", "energy", x);
            break;
        case Energy_Range_Calculator_t::FromRtoE:
            return table.getVal("energy", "range", x);
            break;
    }
}

float Energy_Range_Calculator_t::operator()(float& x)
{
    switch(dir)
    {
        case Energy_Range_Calculator_t::FromEtoR:
            return table.getVal("range", "energy", x);
            break;
        case Energy_Range_Calculator_t::FromRtoE:
            return table.getVal("energy", "range", x);
            break;
    }
}

//std::vector<float> r80 = { 30.710, 33.260, 35.820, 38.360, 40.920, 43.480, 45.980,
//        48.510, 51.030, 53.560, 56.100, 58.630, 61.150, 63.660, 66.190, 68.710, 71.220,
//        73.760, 76.250, 78.770, 81.270, 83.770, 86.280, 88.820, 91.370, 93.950, 96.540,
//        99.100, 101.63, 104.13, 106.61, 109.04, 111.49, 113.94, 116.49, 119.07, 121.69,
//        124.26, 126.71, 129.13, 131.57, 134.11, 136.70, 139.30, 141.89, 144.36, 146.84,
//        149.31, 151.77, 154.31, 156.84, 159.35, 161.88, 164.37, 166.86, 169.37, 171.85,
//        174.35, 176.89, 179.42, 182.02, 184.57, 187.10, 189.56, 192.01, 194.42, 196.87,
//        199.37, 201.90, 204.52, 207.06, 209.58, 212.11, 214.60, 217.07, 219.53, 222.06,
//        224.51, 227.08, 229.63, 232.22, 234.83, 237.41, 240.01, 242.54, 245.05, 247.47,
//        249.88, 252.34, 254.81, 257.24, 259.73, 262.26, 264.75, 267.28, 269.84, 272.38,
//        274.87, 277.32, 279.77, 282.27, 284.78, 287.39, 289.98, 292.57, 295.11, 297.72,
//        300.24, 302.70, 305.16, 307.69, 310.16, 312.68, 315.19, 317.72, 320.23, 322.75,
//        325.23, 327.74, 330.26 };
//std::vector<float> energyJanni = { 59.7690, 62.4880, 65.1150, 67.6600, 70.1300, 72.5360,
//        74.8470, 77.1120, 79.3290, 81.5050, 83.6400, 85.7270, 87.7750, 89.7840, 91.7600,
//        93.7060, 95.6250, 97.5190, 99.3710, 101.201, 102.992, 104.765, 106.527, 108.295,
//        110.042, 111.782, 113.514, 115.209, 116.873, 118.484, 120.069, 121.622, 123.162,
//        124.699, 126.267, 127.861, 129.451, 130.991, 132.461, 133.900, 135.334, 136.811,
//        138.310, 139.814, 141.287, 142.701, 144.096, 145.470, 146.845, 148.235, 149.624,
//        150.998, 152.361, 153.711, 155.041, 156.374, 157.680, 158.993, 160.322, 161.635,
//        162.976, 164.287, 165.584, 166.844, 168.081, 169.294, 170.512, 171.760, 173.026,
//        174.301, 175.561, 176.792, 178.020, 179.226, 180.415, 181.607, 182.799, 183.982,
//        185.191, 186.403, 187.628, 188.844, 190.052, 191.254, 192.421, 193.568, 194.685,
//        195.798, 196.913, 198.040, 199.146, 200.258, 201.388, 202.519, 203.637, 204.774,
//        205.889, 206.982, 208.073, 209.144, 210.230, 211.333, 212.454, 213.577, 214.699,
//        215.795, 216.890, 217.962, 219.026, 220.075, 221.131, 222.177, 223.237, 224.291,
//        225.345, 226.389, 227.418, 228.460, 229.484, 230.527 };

// DataTable_t Janni2ICRU[3][77] =
// {
//     // Energies MeV
//     { 0.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,12.5,15.0,17.5,20.0,22.5,25.0,27.5,30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 350, 400 },
//     // Janni82 range mm water
//     { 0.0, 0.3059, 0.3676, 0.4342, 0.5058, 0.5822, 0.6634, 0.7493, 0.8399, 0.9351, 1.0348, 1.1391, 1.2479, 1.8575, 2.5739, 3.3939, 4.3147, 5.3336, 6.4487, 7.6579, 8.9593, 11.833, 15.057, 18.620, 22.515, 26.730, 31.259, 36.093, 41.226, 46.650, 52.360, 58.348, 64.610, 71.138, 77.929, 84.960, 92.263, 99.802, 107.57, 115.64, 123.84, 132.33, 141.08, 150.03, 159.19, 168.49, 178.20, 187.92, 197.92, 208.07, 218.44, 228.99, 239.74, 250.72, 261.85, 273.18, 284.65, 296.34, 308.07, 320.13, 332.21, 344.54, 357.00, 369.73, 382.57, 395.41, 408.58, 421.82, 435.30, 448.86, 462.64, 476.44, 490.41, 504.42, 518.71, 668.12, 829.02 },
//     // ICRU49 range mm water
//     { 0.0, 0.3007, 0.3614, 0.4268, 0.4972, 0.5734, 0.6522, 0.7368, 0.8259, 0.9196, 1.0179, 1.1206, 1.2275, 1.8285, 2.5344, 3.3440, 4.2530, 5.2580, 6.3590, 7.5540, 8.8390, 11.681, 14.868, 18.382, 22.237, 26.400, 30.887, 35.670, 40.743, 46.115, 51.767, 57.689, 63.890, 70.351, 77.080, 84.083, 91.267, 98.660, 106.42, 114.45, 122.51, 131.03, 139.54, 148.44, 157.50, 166.80, 176.29, 185.90, 195.94, 205.95, 216.27, 226.75, 237.24, 248.16, 259.29, 270.48, 281.92, 293.33, 305.10, 317.02, 329.11, 341.33, 353.71, 364.94, 378.98, 391.74, 404.68, 417.98, 431.43, 444.71, 458.40, 471.86, 485.87, 499.81, 513.93, 662.07, 821.68 }
// };

// std::vector< std::vector<double> > R80AstroidTopasb8[3][27] =
// {
//     // Energy [MeV]
//     { 91.015, 95.489, 101.50, 109.36, 117.30, 122.69, 127.81, 134.27, 139.28, 146.68, 152.35, 156.37, 161.56, 166.79, 171.42, 176.88, 181.97, 186.06, 190.47, 195.37, 199.13, 204.11, 208.37, 212.67, 216.58, 221.15, 223.58 },
//     // R80 astroid [mm]
//     { 65.420, 71.310, 79.590, 91.020, 103.06, 111.54, 119.84, 130.77, 139.50, 152.79, 163.30, 170.89, 180.91, 191.30, 200.60, 211.82, 222.48, 231.22, 240.77, 251.57, 259.96, 271.24, 281.02, 291.04, 300.26, 311.14, 317.04 },
//     // R80 TOPASb8 [mm]
//     { 65.040, 70.850, 78.980, 90.170, 102.08, 110.52, 118.77, 129.54, 138.14, 151.25, 161.62, 169.13, 179.04, 189.24, 198.44, 209.53, 220.06, 228.67, 238.09, 248.74, 257.02, 268.16, 277.82, 287.70, 296.79, 307.55, 313.34 }
// };

// std::vector< std::vector<double> > R80R90EnergyMGHR1[3][129] =
// {
//     // R90 mm water
//     { 65.080, 67.860, 70.930, 73.510, 76.320, 79.160, 82.000, 84.860, 87.720, 90.550, 93.990, 96.830, 99.690, 102.53, 105.37, 108.19, 110.99, 113.72, 116.48, 119.25, 121.99, 124.69, 127.40, 130.11, 132.80, 135.48, 138.81, 141.51, 144.20, 146.85, 149.43, 152.06, 154.71, 157.31, 159.94, 162.52, 165.05, 167.58, 170.10, 172.62, 175.13, 177.64, 180.10, 182.50, 185.59, 188.03, 190.45, 192.82, 195.12, 197.45, 199.73, 202.01, 204.33, 206.53, 208.74, 210.94, 213.09, 215.24, 217.41, 219.54, 221.57, 223.64, 225.72, 228.30, 230.31, 232.28, 234.18, 236.07, 237.95, 239.84, 241.74, 243.58, 245.36, 247.14, 248.90, 250.64, 252.36, 254.07, 255.73, 257.41, 259.01, 260.59, 262.18, 263.73, 265.79, 267.31, 268.82, 270.29, 271.73, 273.14, 274.55, 275.95, 277.34, 278.72, 280.07, 281.34, 282.62, 283.90, 285.15, 286.41, 287.67, 288.88, 290.09, 291.66, 292.78, 293.89, 295.01, 296.10, 297.17, 298.23, 299.31, 300.35, 301.38, 302.42, 303.40, 304.36, 305.32, 306.26, 307.18, 308.09, 309.29, 310.18, 311.07, 311.95, 312.82, 313.66, 314.49, 315.27, 316.08 },
//     // R80 mm water
//     { 65.420, 68.220, 71.320, 73.900, 76.730, 79.590, 82.440, 85.320, 88.180, 91.030, 94.480, 97.340, 100.21, 103.06, 105.92, 108.75, 111.56, 114.30, 117.07, 119.86, 122.61, 125.32, 128.04, 130.77, 133.47, 136.15, 139.51, 142.21, 144.91, 147.57, 150.16, 152.80, 155.45, 158.07, 160.70, 163.29, 165.83, 168.37, 170.90, 173.43, 175.94, 178.46, 180.92, 183.34, 186.44, 188.88, 191.31, 193.68, 195.99, 198.32, 200.61, 202.89, 205.21, 207.43, 209.64, 211.82, 213.98, 216.15, 218.32, 220.45, 222.49, 224.57, 226.63, 229.22, 231.23, 233.20, 235.11, 237.00, 238.89, 240.78, 242.68, 244.52, 246.31, 248.09, 249.84, 251.58, 253.31, 255.02, 256.68, 258.36, 259.96, 261.54, 263.14, 264.68, 266.74, 268.27, 269.78, 271.24, 272.67, 274.10, 275.51, 276.91, 278.30, 279.67, 281.02, 282.29, 283.58, 284.87, 286.12, 287.38, 288.63, 289.84, 291.05, 292.62, 293.74, 294.85, 295.97, 297.05, 298.13, 299.18, 300.27, 301.31, 302.34, 303.38, 304.36, 305.32, 306.27, 307.21, 308.14, 309.05, 310.25, 311.15, 312.04, 312.92, 313.78, 314.62, 315.45, 316.24, 317.04 },
//     // energy
//     { 91.015, 93.164, 95.489, 97.421, 99.481, 101.50, 103.48, 105.47, 107.43, 109.36, 111.67, 113.56, 115.45, 117.30, 119.13, 120.93, 122.69, 124.41, 126.11, 127.81, 129.46, 131.08, 132.68, 134.27, 135.83, 137.37, 139.28, 140.80, 142.31, 143.78, 145.22, 146.68, 148.12, 149.53, 150.96, 152.35, 153.70, 155.04, 156.37, 157.69, 159.00, 160.30, 161.56, 162.79, 164.35, 165.58, 166.79, 167.98, 169.13, 170.29, 171.42, 172.54, 173.67, 174.75, 175.82, 176.88, 177.92, 178.96, 179.99, 181.01, 181.97, 182.95, 183.91, 185.12, 186.06, 186.98, 187.86, 188.74, 189.61, 190.47, 191.33, 192.17, 192.98, 193.79, 194.58, 195.37, 196.15, 196.92, 197.66, 198.42, 199.13, 199.83, 200.55, 201.23, 202.14, 202.81, 203.47, 204.11, 204.74, 205.37, 205.98, 206.59, 207.19, 207.79, 208.37, 208.93, 209.48, 210.04, 210.57, 211.11, 211.64, 212.15, 212.67, 213.34, 213.81, 214.29, 214.77, 215.23, 215.68, 216.13, 216.58, 217.02, 217.45, 217.89, 218.30, 218.71, 219.11, 219.50, 219.89, 220.27, 220.78, 221.15, 221.52, 221.88, 222.24, 222.58, 222.92, 223.25, 223.58 }
// };

// std::vector< std::vector<double> > R80R90EnergyTopasA5[3][120] =
// {
//     // R90 mm water
//     { 30.550, 33.090, 35.640, 38.170, 40.720, 43.260, 45.760, 48.280, 50.780, 53.310, 55.830, 58.360, 60.870, 63.380, 65.900, 68.390, 70.910, 73.430, 75.910, 78.420, 80.910, 83.390, 85.900, 88.430, 90.990, 93.550, 96.130, 98.680, 101.20, 103.67, 106.16, 108.55, 111.02, 113.46, 115.99, 118.58, 121.17, 123.73, 126.17, 128.59, 131.00, 133.53, 136.12, 138.70, 141.29, 143.75, 146.23, 148.67, 151.13, 153.63, 156.18, 158.67, 161.21, 163.71, 166.15, 168.68, 171.13, 173.60, 176.15, 178.67, 181.25, 183.80, 186.29, 188.76, 191.17, 193.64, 196.06, 198.55, 201.06, 203.68, 206.19, 208.73, 211.26, 213.73, 216.18, 218.64, 221.17, 223.60, 226.14, 228.72, 231.28, 233.91, 236.44, 239.07, 241.60, 244.06, 246.48, 248.90, 251.40, 253.85, 256.24, 258.71, 261.25, 263.76, 266.26, 268.82, 271.36, 273.81, 276.27, 278.74, 281.23, 283.73, 286.34, 288.96, 291.48, 294.03, 296.68, 299.18, 301.63, 304.09, 306.59, 309.11, 311.59, 314.10, 316.66, 319.15, 321.66, 324.11, 326.65, 329.15 },
//     // R80 mm water
//     { 30.710, 33.260, 35.820, 38.360, 40.920, 43.480, 45.980, 48.510, 51.030, 53.560, 56.100, 58.630, 61.150, 63.660, 66.190, 68.710, 71.220, 73.760, 76.250, 78.770, 81.270, 83.770, 86.280, 88.820, 91.370, 93.950, 96.540, 99.100, 101.63, 104.13, 106.61, 109.04, 111.49, 113.94, 116.49, 119.07, 121.69, 124.26, 126.71, 129.13, 131.57, 134.11, 136.70, 139.30, 141.89, 144.36, 146.84, 149.31, 151.77, 154.31, 156.84, 159.35, 161.88, 164.37, 166.86, 169.37, 171.85, 174.35, 176.89, 179.42, 182.02, 184.57, 187.10, 189.56, 192.01, 194.42, 196.87, 199.37, 201.90, 204.52, 207.06, 209.58, 212.11, 214.60, 217.07, 219.53, 222.06, 224.51, 227.08, 229.63, 232.22, 234.83, 237.41, 240.01, 242.54, 245.05, 247.47, 249.88, 252.34, 254.81, 257.24, 259.73, 262.26, 264.75, 267.28, 269.84, 272.38, 274.87, 277.32, 279.77, 282.27, 284.78, 287.39, 289.98, 292.57, 295.11, 297.72, 300.24, 302.70, 305.16, 307.69, 310.16, 312.68, 315.19, 317.72, 320.23, 322.75, 325.23, 327.74, 330.26 },
//     // energy
//     { 59.7690, 62.4880, 65.1150, 67.6600, 70.1300, 72.5360, 74.8470, 77.1120, 79.3290, 81.5050, 83.6400, 85.7270, 87.7750, 89.7840, 91.7600, 93.7060, 95.6250, 97.5190, 99.3710, 101.201, 102.992, 104.765, 106.527, 108.295, 110.042, 111.782, 113.514, 115.209, 116.873, 118.484, 120.069, 121.622, 123.162, 124.699, 126.267, 127.861, 129.451, 130.991, 132.461, 133.900, 135.334, 136.811, 138.310, 139.814, 141.287, 142.701, 144.096, 145.470, 146.845, 148.235, 149.624, 150.998, 152.361, 153.711, 155.041, 156.374, 157.680, 158.993, 160.322, 161.635, 162.976, 164.287, 165.584, 166.844, 168.081, 169.294, 170.512, 171.760, 173.026, 174.301, 175.561, 176.792, 178.020, 179.226, 180.415, 181.607, 182.799, 183.982, 185.191, 186.403, 187.628, 188.844, 190.052, 191.254, 192.421, 193.568, 194.685, 195.798, 196.913, 198.040, 199.146, 200.258, 201.388, 202.519, 203.637, 204.774, 205.889, 206.982, 208.073, 209.144, 210.230, 211.333, 212.454, 213.577, 214.699, 215.795, 216.890, 217.962, 219.026, 220.075, 221.131, 222.177, 223.237, 224.291, 225.345, 226.389, 227.418, 228.460, 229.484, 230.527 }
// };
