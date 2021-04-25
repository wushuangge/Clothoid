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

    double getCoefficient(double startCurvature, double endCurvature, double length, int sign, Coordinates c1, Coordinates c2){
        double ASQ = fabs(length / (endCurvature - startCurvature));
        double La = ASQ * startCurvature;
        double Le = ASQ * endCurvature;

        double tempx0 = integral(cosFunction, 0, La, ASQ);
        double tempy0 = integral(sinFunction, 0, La, ASQ);

        double tempx1 = integral(cosFunction, 0, Le, ASQ);
        double tempy1 = integral(sinFunction, 0, Le, ASQ);

        double a1 = atan2(c2.y - c1.y, c2.x - c1.x);
        double a2 = atan2(sign * tempy1 - sign * tempy0, tempx1 - tempx0);
        //double a3 = atan2(sign * tempy0 - sign * tempy1, tempx0 - tempx1);
        double angle = a1 - a2;

        double X00  = (tempx0 * cos(angle) - (sign * tempy0) * sin(angle));
        double Y00  = (tempx0 * sin(angle) + (sign * tempy0) * cos(angle));

        double X01  = (tempx1 * cos(angle) - (sign * tempy1) * sin(angle));
        double Y01  = (tempx1 * sin(angle) + (sign * tempy1) * cos(angle));

        double X0000 = X01 - X00;
        double Y0000 = Y01 - Y00;

        double X0001 = c2.x - c1.x;
        double Y0001 = c2.y - c1.y;

        double coefficient = (X0001 / X0000 +  Y0001 / Y0000) / 2;
        return coefficient;
    }

    Coordinates calcOriginCoordinates(int originX, int originY, int dx, int dy){
        Coordinates coordinates;
        coordinates.x = originX + dx * 32;
        coordinates.y = originY + dy * 32;
        std::cout<<"coordinates: " << std::setprecision(15) <<coordinates.x<<";"<<std::setprecision(15) <<coordinates.y << std::endl;
        return coordinates;
    }

    int getSign(double startCurvature, double endCurvature, double tangentArcStart, double tangentArcEnd){
        int sign = 0;
        double angleStart = 90.0f - ((double)tangentArcStart * 360.0f / 65535);
        while (angleStart < 0)
            angleStart = angleStart + 360;
        double angleEnd = 90.0f - ((double)tangentArcEnd * 360.0f / 65535);
        while (angleEnd < 0)
            angleEnd = angleEnd + 360;
        std::cout<< "tangentArcStart: "<<angleStart<<std::endl;
        std::cout<< "tangentArcEnd: "<<angleEnd<<std::endl;

        if (startCurvature < endCurvature && angleStart > angleEnd)
            sign = -1;
        else if (startCurvature < endCurvature && angleStart < angleEnd)
            sign = 1;
        else if (startCurvature > endCurvature && (angleStart > angleEnd))
            sign = -1;
        else if (startCurvature > endCurvature && (angleStart < angleEnd))
            sign = 1;
        std::cout<< "sign: "<<sign<<std::endl;
        return sign;
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
