//
// Created by wusg on 2021/4/21.
//

#ifndef MXCLOTHOID_MXCLOTHOID_H
#define MXCLOTHOID_MXCLOTHOID_H

#include <cmath>
#include <iomanip>
#include <vector>
#include "Clothoid.h"

#define PI 3.14159265358979323846264338328
#define ARRAY_SIZE(x) (sizeof(x)/sizeof((x)[0]))
#define ANGLE_DIR(a, b) ((a >= 0.0f) ? 0 : (b == 0 ? 0 : 1))    //0: right, 1, left
#define ANGLE_SIGN(a) ((a == 0) ? -1 : 1)

namespace MXClothoid{

    struct Coordinates{
        double x;
        double y;
    };

    // A=Ca  B=Cb-Ca/(2*L) L=回旋曲线总长度 l=任意一点到曲线起点的长度
    double ClothoidX(double A, double B,double l){
        double x = l - pow(A,2)*pow(l,3)/6 - A*B*pow(l,4)/4 - (pow(B,2)/10-pow(A,4)/120)*pow(l,5) +
                   pow(A,3)*B*pow(l,6)/36 + pow(A,2)*pow(B,2)*pow(l,7)/28 + A*pow(B,3)*pow(l,8)/48 + pow(B,4)*pow(l,9)/216;
        return x;
    }

    // A=Ca  B=Cb-Ca/(2*L) L=回旋曲线总长度 l=任意一点到曲线起点的长度
    double ClothoidY(double A, double B,double l){
        double y = A*pow(l,2)/2 + B*pow(l,3)/3 - pow(A,3)*pow(l,4)/24 -
                   pow(A,2)*B*pow(l,5)/10 - (A*pow(B,2)/12 - pow(A,5)/720)*pow(l,6) -
                   (pow(B,3)/42 - pow(A,4)*B/168)*pow(l,7) +pow(A,3)*pow(B,2)*pow(l,8)/96 + pow(A,2)*pow(B,3)*pow(l,9)/108;
        return y;
    }

    double cosFunction(double x, double asquare)
    {
        double y;
        y = cos((x * x)/(2 * asquare));
        return y;
    }

    double sinFunction(double x, double asquare)
    {
        double y;
        y = sin((x * x)/(2 * asquare));
        return y;
    }

    double integral(double(*f)(double x, double asquare), double a, double b, double asquare, int n = 400) {
        double step = (b - a) / n;  // width of each small rectangle
        double area = 0.0;  // signed area
        for (int i = 0; i < n; i ++) {
            area += f(a + (i + 0.5) * step, asquare) * step; // sum up each small rectangle
        }
        return area;
    }

    template <typename T> int SIGN(T val) {
        return (T(0) < val) - (val < T(0));
    }

    double getcurvature(int curvature){
        int C = curvature - 511;
        double temp = 0.0;
        if (C == 0)
            temp = 0.0;
        if (1 <= abs(C) && abs(C) < 64)
            temp = (double)512 * (double)(C + SIGN(C) * 19) / (double)1038539554;
        if (64 <= abs(C) && abs(C) < 128)
            temp = (double)2048 * (double)(C - SIGN(C) * 44) / (double)988894254;
        if (128 <= abs(C) && abs(C) < 192)
            temp = (double)8192 * (double)(C - SIGN(C) * 108) / (double)932923357;
        if (192 <= abs(C) && abs(C) < 256)
            temp = (double)32768 * (double)(C - SIGN(C) * 172) / (double)879961330;
        if (256 <= abs(C) && abs(C) < 320)
            temp = (double)131072 * (double)(C - SIGN(C) * 236) / (double)829979357;
        if (320 <= abs(C) && abs(C) < 384)
            temp = (double)524288 * (double)(C - SIGN(C) * 300) / (double)782847732;
        if (384 <= abs(C) && abs(C) < 448)
            temp = (double)2097152 * (double)(C - SIGN(C) * 364) / (double)738390903;
        if (448 <= abs(C) && abs(C) <= 511)
            temp = (double)8388608 * (double)(C - SIGN(C) * 428) / (double)696254464;
        return temp;
    }

    Coordinates calcOriginCoordinates(int originX, int originY, int dx, int dy){
        Coordinates coordinates;
        coordinates.x = originX + dx * 32;
        coordinates.y = originY + dy * 32;
        std::cout<<"coordinates: " << std::setprecision(15) <<coordinates.x<<";"<<std::setprecision(15) <<coordinates.y << std::endl;
        return coordinates;
    }

    int convertCurvatureToProtocol(float curvature) {
        static const float sgCurvatureConvertTable[] = {
                0.9329005,0.7415425,0.5174103,0.3524464,0.2462423,0.1785520,0.1341612,0.1039609,
                0.0826743,0.0671898,0.0556131,0.0467509,0.0398266,0.0343195,0.0298710,0.0262284,
                0.0232093,0.0206800,0.0185407,0.0167154,0.0151458,0.0137865,0.0126017,0.0115628,
                0.0106469,0.0098354,0.0091130,0.0084671,0.0078875,0.0073652,0.0068931,0.0064649,
                0.0060753,0.0057198,0.0053945,0.0050962,0.0048220,0.0045693,0.0043359,0.0041200,
                0.0039197,0.0037337,0.0035606,0.0033993,0.0032487,0.0031078,0.0029759,0.0028523,
                0.0027362,0.0026270,0.0025242,0.0024273,0.0023359,0.0022496,0.0021679,0.0020907,
                0.0020175,0.0019480,0.0018821,0.0018195,0.0017599,0.0017032,0.0016492,0.0015978,
                0.0015487,0.0015018,0.0014570,0.0014142,0.0013733,0.0013341,0.0012966,0.0012606,
                0.0012261,0.0011930,0.0011612,0.0011307,0.0011014,0.0010732,0.0010460,0.0010199,
                0.0009947,0.0009705,0.0009471,0.0009246,0.0009029,0.0008819,0.0008616,0.0008420,
                0.0008231,0.0008049,0.0007872,0.0007701,0.0007536,0.0007375,0.0007220,0.0007070,
                0.0006924,0.0006783,0.0006646,0.0006513,0.0006384,0.0006259,0.0006138,0.0006020,
                0.0005905,0.0005794,0.0005686,0.0005580,0.0005478,0.0005378,0.0005282,0.0005188,
                0.0005096,0.0005006,0.0004919,0.0004834,0.0004752,0.0004671,0.0004593,0.0004516,
                0.0004442,0.0004369,0.0004298,0.0004229,0.0004161,0.0004095,0.0004030,0.0003967,
                0.0003906,0.0003846,0.0003787,0.0003730,0.0003674,0.0003619,0.0003566,0.0003513,
                0.0003462,0.0003412,0.0003363,0.0003315,0.0003268,0.0003222,0.0003177,0.0003132,
                0.0003089,0.0003047,0.0003006,0.0002965,0.0002925,0.0002886,0.0002848,0.0002811,
                0.0002774,0.0002738,0.0002703,0.0002668,0.0002634,0.0002601,0.0002568,0.0002536,
                0.0002504,0.0002473,0.0002443,0.0002413,0.0002384,0.0002355,0.0002327,0.0002299,
                0.0002272,0.0002245,0.0002219,0.0002193,0.0002168,0.0002143,0.0002118,0.0002094,
                0.0002071,0.0002048,0.0002025,0.0002002,0.0001980,0.0001958,0.0001937,0.0001916,
                0.0001895,0.0001875,0.0001855,0.0001835,0.0001816,0.0001797,0.0001778,0.0001759,
                0.0001741,0.0001723,0.0001705,0.0001688,0.0001671,0.0001654,0.0001637,0.0001621,
                0.0001605,0.0001589,0.0001573,0.0001558,0.0001543,0.0001528,0.0001513,0.0001498,
                0.0001484,0.0001470,0.0001456,0.0001442,0.0001429,0.0001416,0.0001403,0.0001390,
                0.0001377,0.0001364,0.0001351,0.0001339,0.0001327,0.0001315,0.0001303,0.0001292,
                0.0001280,0.0001269,0.0001258,0.0001247,0.0001236,0.0001225,0.0001214,0.0001204,
                0.0001194,0.0001183,0.0001173,0.0001163,0.0001153,0.0001144,0.0001134,0.0001125,
                0.0001116,0.0001107,0.0001097,0.0001088,0.0001080,0.0001071,0.0001062,0.0001053,
                0.0001045,0.0001036,0.0001028,0.0001020,0.0001012,0.0001004,0.0000500,0.0000000,
        };

        int l = -1;
        int r = ARRAY_SIZE(sgCurvatureConvertTable);


        while (l < r - 1) {
            int m = (l + r) / 2;
            if (curvature < sgCurvatureConvertTable[m]) {
                l = m;
            } else {
                r = m;
            }
        }
        return (r >= 0 && r < ARRAY_SIZE(sgCurvatureConvertTable)) ? r : 255;
    }

    double calcBranchAngle(int32_t angle_A, int32_t angle_B){
        float branch_angle = (angle_B - angle_A) * 360.0f / 65535;
        if (branch_angle > 180.0f) {
            return branch_angle - 360.0f;
        } else if (branch_angle < -180.0f) {
            return branch_angle + 360.0f;
        }
        return branch_angle;
    }

    double calcEndTangentArc(double startCurvature, double endCurvature, double length, int sign, double tangentArc){
        double asquare = fabs(length / (endCurvature - startCurvature));
        double La = asquare * startCurvature;
        double Le = asquare * endCurvature;

        double beta1 = sign * pow(La, 2)/(2 * asquare) * 180.0f / PI;     //p 点切线方位角
        while (beta1 < 0)
            beta1 = beta1 + 360.0f;

        double beta2 = sign * pow(Le, 2)/(2 * asquare) * 180.0f / PI;     //p 点切线方位角
        while (beta2 < 0)
            beta2 = beta2 + 360.0f;

        std::cout<< "beta_diff: "<<beta2 - beta1<<std::endl;
        double ang = (90.0f - ((double)tangentArc * 360.0f / 65535));
        while (ang < 0)
            ang = ang + 360.0f;

        double angle = ang - beta1 + beta2;
        while (angle < 0)
            angle = angle + 360.0f;
        return angle;
    }

    double calcStartAngle(double startCurvature, double endCurvature, double length, int sign, Coordinates c1, Coordinates c2){
        double asquare = fabs(length / (endCurvature - startCurvature));
        double La = asquare * startCurvature;
        double Le = asquare * endCurvature;
        double tempx = MXClothoid::integral(MXClothoid::cosFunction, La, Le, asquare);
        double tempy = MXClothoid::integral(MXClothoid::sinFunction, La, Le, asquare);

        double a1 = atan2(c2.y - c1.y, c2.x - c1.x);
        double a2= atan2(sign * tempy, tempx);
        double angle = a1 - a2;     //a1 - a2 表示 独立坐标按线路旋转的角度
        while (angle < 0)
            angle = angle + 2 * PI;
        std::cout<<"angle: "<<(angle*180/PI)<<std::endl;
        return angle;
    }

    double calcEndAngle(double startCurvature, double endCurvature, double length, int sign, Coordinates c1, Coordinates c2){
        double asquare = fabs(length / (endCurvature - startCurvature));
        double La = asquare * startCurvature;
        double Le = asquare * endCurvature;
        double tempx = MXClothoid::integral(MXClothoid::cosFunction, Le, La, asquare);
        double tempy = MXClothoid::integral(MXClothoid::sinFunction, Le, La, asquare);

        double tempx1 = MXClothoid::integral(MXClothoid::cosFunction, 0, La, asquare);
        double tempy1 = MXClothoid::integral(MXClothoid::sinFunction, 0, La, asquare);

        double tempx2 = MXClothoid::integral(MXClothoid::cosFunction, 0, Le, asquare);
        double tempy2 = MXClothoid::integral(MXClothoid::sinFunction, 0, Le, asquare);

        double a1 = atan2(c2.y - c1.y, c2.x - c1.x);
        double a2 = atan2(sign * tempy, tempx);
        double a3 = atan2(sign * tempy2 - sign * tempy1, tempx2 - tempx1);
        double angle = a1 - a3;     //a1 - a2 表示 独立坐标按线路旋转的角度
        while (angle < 0)
            angle = angle + 2 * PI;
        std::cout<<"angle: "<<(angle*180/PI)<<std::endl;
        return angle;
    }

    double calcAngleByStartTangentArc(double startCurvature, double endCurvature, double length, int sign, double tangentArc){
        double ang = (90.0f - ((double)tangentArc * 360.0f / 65535)) * PI / 180.0f;
        double asquare = fabs(length / (endCurvature - startCurvature));
        double La = asquare * startCurvature;
        double beta = sign * pow(La, 2)/(2 * asquare);
        double angle = ang - beta;
        while (angle < 0)
            angle = angle + 2 * PI;
        std::cout<<"angle: "<<(angle*180/PI)<<std::endl;
        return angle;
    }

    double calcAngleByEndTangentArc(double startCurvature, double endCurvature, double length, int sign, double tangentArc){
        double ang = (90.0f - ((double)tangentArc * 360.0f / 65535)) * PI / 180.0f + PI;
        double asquare = fabs(length / (endCurvature - startCurvature));
        double Le = asquare * endCurvature;
        double beta = sign * pow(Le, 2)/(2 * asquare);
        while (beta < 0)
            beta = beta + 2 * PI;
        double angle = ang - beta;
        while (angle < 0)
            angle = angle + 2 * PI;
        std::cout<<"angle: "<<(angle*180/PI)<<std::endl;
        return angle;
    }

    std::vector<Coordinates> calcClothoidCoordinates(double startCurvature, double endCurvature, double length, int npts, int sign, double angle){
        std::vector<Coordinates> veccoord;
        if(startCurvature == endCurvature){
            std::cout<<"起点曲率不能等于终点曲率!!!"<<std::endl;
            return veccoord;
        }
        double asquare = fabs(length / (endCurvature - startCurvature));
        double La = asquare * startCurvature;
        double Le = asquare * endCurvature;
        double averageLength =  length / npts;
        for (int i = 0; i <= npts; i++) {
            double a = 0.0;
            double b = 0.0;
            if(La > Le){
                a = Le;
                b = Le + averageLength * i;
            }else{
                a = La;
                b = La + averageLength * i;
            }
            double tempx = integral(cosFunction, a, b, asquare);
            double tempy = integral(sinFunction, a, b, asquare);

            Coordinates coord;
            coord.x  = ((tempx * cos(angle) - (sign * tempy) * sin(angle)));
            coord.y  = ((tempx * sin(angle) + (sign * tempy) * cos(angle)));

            veccoord.push_back(coord);
        }
        return veccoord;
    }
}

#endif //MXCLOTHOID_MXCLOTHOID_H
