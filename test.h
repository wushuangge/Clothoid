//
// Created by wusg on 2021/4/21.
//

#ifndef MXCLOTHOID_TEST_H
#define MXCLOTHOID_TEST_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>
#include "Clothoid.h"
#include "MXClothoid.h"

namespace testClothoid{
    //验证计算开源库
    std::vector<MXClothoid::Coordinates> calcClothoid(int npts, double length, double tangentArcStart, double tangentArcEnd, MXClothoid::Coordinates coord1, MXClothoid::Coordinates coord2){
        std::vector<MXClothoid::Coordinates> vecCoordinates;
        double psd_calc_k = 0.0;
        double psd_calc_dk = 0.0;
        double psd_calc_l = 0.0;

        Clothoid::buildClothoid(
                (double)coord1.x,
                (double)coord1.y,
                (double)tangentArcStart * 360.0f / 65535,
                (double)coord2.x,
                (double)coord2.y,
                (double)tangentArcEnd * 360.0f / 65535,
                psd_calc_k, psd_calc_dk, psd_calc_l);

        double averageLength = length / npts;
        for (int i = 1; i <= npts; ++i) {
            double X = 0.0;
            double Y = 0.0;
            Clothoid::pointsOnClothoid(
                    (double)coord1.x,
                    (double)coord1.y,
                    (double)tangentArcStart * 360.0f / 65535,
                    (double)psd_calc_k,
                    (double)psd_calc_dk,
                    (double)averageLength * i,
                    X,Y);
            MXClothoid::Coordinates coordtemp;
            coordtemp.x = X;
            coordtemp.y = Y;
            vecCoordinates.push_back(coordtemp);
        }
        return vecCoordinates;
    }

    ///计算地图经纬度比例
    double getCoefficient(double startCurvature, double endCurvature, double length, int sign, MXClothoid::Coordinates c1, MXClothoid::Coordinates c2){
        double ASQ = fabs(length / (endCurvature - startCurvature));
        double La = ASQ * startCurvature;
        double Le = ASQ * endCurvature;

        double tempx0 = MXClothoid::integral(MXClothoid::cosFunction, 0, La, ASQ);
        double tempy0 = MXClothoid::integral(MXClothoid::sinFunction, 0, La, ASQ);

        double tempx1 = MXClothoid::integral(MXClothoid::cosFunction, 0, Le, ASQ);
        double tempy1 = MXClothoid::integral(MXClothoid::sinFunction, 0, Le, ASQ);

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

    //通过起始点坐标进行计算
    void calcCoordinates1(double startCurvature, double endCurvature, double tangentArcStart,
                          double tangentArcEnd, MXClothoid::Coordinates c1, MXClothoid::Coordinates c2,
                          double length, int npts){
        //计算左旋还是右旋
        int branchAngleDirection = 0;
        if (startCurvature < endCurvature){
            //第一缓和曲线
            //计算角度差
            double branch_angle = MXClothoid::calcBranchAngle(tangentArcStart, tangentArcEnd);
            int tmp = MXClothoid::convertCurvatureToProtocol(startCurvature);
            branchAngleDirection = ANGLE_DIR(branch_angle, tmp);  //0: right, 1, left
        }else{
            //第二缓和曲线
            //计算角度差
            double branch_angle = MXClothoid::calcBranchAngle(tangentArcEnd, tangentArcStart);
            int tmp = MXClothoid::convertCurvatureToProtocol(endCurvature);
            branchAngleDirection = ANGLE_DIR(branch_angle, tmp);  //0: right, 1, left
        }
        int sign = ANGLE_SIGN(branchAngleDirection);
        //计算地图缩放比例
        double coefficient = getCoefficient(startCurvature, endCurvature, length, sign, c1, c2);

        std::vector<MXClothoid::Coordinates> veccoord;
        double angle = 0.0;
        std::ofstream outFile;
        outFile.open("../data.csv", std::ios::app|std::ios::out);
        if (startCurvature < endCurvature){
            //计算独立坐标相对偏转角度,用于坐标转换
            angle = MXClothoid::calcStartAngle(startCurvature, endCurvature, length, sign, c1, c2);
            veccoord = MXClothoid::calcClothoidCoordinates(startCurvature, endCurvature, length, npts, sign, angle);
            for (int i = 0; i < veccoord.size(); i++) {
                double x = c1.x + veccoord[i].x * coefficient;
                double y = c1.y + veccoord[i].y * coefficient;
                outFile << std::setprecision(15)<< x << ',' << std::setprecision(15)<< y << std::endl;
                std::cout<<i <<","<< std::setprecision(15) << x <<","<<std::setprecision(15) << y << "," << i << std::endl;
            }
        }
        else{
            //计算独立坐标相对偏转角度,用于坐标转换
            angle = MXClothoid::calcEndAngle(startCurvature, endCurvature, length, sign, c1, c2);
            veccoord = MXClothoid::calcClothoidCoordinates(startCurvature, endCurvature, length, npts, sign, angle);
            for (int i = 0; i < veccoord.size(); i++) {
                double x = c2.x + veccoord[i].x * coefficient;
                double y = c2.y + veccoord[i].y * coefficient;
                outFile << std::setprecision(15)<< x << ',' << std::setprecision(15)<< y << std::endl;
                std::cout<<i <<","<< std::setprecision(15) << x <<","<<std::setprecision(15) << y << "," << i << std::endl;
            }
        }
        outFile.close();
    }

    ///通过切线方位角计算
    void calcCoordinates2(double startCurvature, double endCurvature, double tangentArcStart,
                          double tangentArcEnd, MXClothoid::Coordinates c1, MXClothoid::Coordinates c2,
                          double length, int npts){
        //计算左旋还是右旋
        int branchAngleDirection = 0;
        if (startCurvature < endCurvature){
            //第一缓和曲线
            //计算角度差
            double branch_angle = MXClothoid::calcBranchAngle(tangentArcStart, tangentArcEnd);
            int tmp = MXClothoid::convertCurvatureToProtocol(startCurvature);
            branchAngleDirection = ANGLE_DIR(branch_angle, tmp);  //0: right, 1, left
        }else{
            //第二缓和曲线
            //计算角度差
            double branch_angle = MXClothoid::calcBranchAngle(tangentArcEnd, tangentArcStart);
            int tmp = MXClothoid::convertCurvatureToProtocol(endCurvature);
            branchAngleDirection = ANGLE_DIR(branch_angle, tmp);  //0: right, 1, left
        }
        int sign = ANGLE_SIGN(branchAngleDirection);
        //计算地图缩放比例
        double coefficient = getCoefficient(startCurvature, endCurvature, length, sign, c1, c2);

        std::vector<MXClothoid::Coordinates> veccoord;
        std::ofstream outFile;
        outFile.open("../data.csv", std::ios::app|std::ios::out);
        if (startCurvature < endCurvature){
            //计算独立坐标相对偏转角度,用于坐标转换
            double startAngle = (90.0f - ((double)tangentArcStart * 360.0f / 65535)) * PI / 180.0f;
            double angle = MXClothoid::calcAngleByStartTangentArc(startCurvature, endCurvature, length, sign, startAngle);
            veccoord = MXClothoid::calcClothoidCoordinates(startCurvature, endCurvature, length, npts, sign, angle);
            for (int i = 0; i < veccoord.size(); i++) {
                double x = c1.x + veccoord[i].x * coefficient;
                double y = c1.y + veccoord[i].y * coefficient;
                outFile << std::setprecision(15)<< x << ',' << std::setprecision(15)<< y << std::endl;
                std::cout<<i <<","<< std::setprecision(15) << x <<","<<std::setprecision(15) << y << "," << i << std::endl;
            }
        }
        else{
            //计算独立坐标相对偏转角度,用于坐标转换
            double endAngle = (90.0f - ((double)tangentArcEnd * 360.0f / 65535)) * PI / 180.0f;
            double angle = MXClothoid::calcAngleByEndTangentArc(startCurvature, endCurvature, length, sign, endAngle);
            veccoord = MXClothoid::calcClothoidCoordinates(startCurvature, endCurvature, length, npts, sign, angle);
            for (int i = 0; i < veccoord.size(); i++) {
                double x = c2.x + veccoord[i].x * coefficient;
                double y = c2.y + veccoord[i].y * coefficient;
                outFile << std::setprecision(15)<< x << ',' << std::setprecision(15)<< y << std::endl;
                std::cout<<i <<","<< std::setprecision(15) << x <<","<<std::setprecision(15) << y << "," << i << std::endl;
            }
        }
        outFile.close();
    }

    //第二缓和曲线也当做第一缓和曲线计算
    void calcCoordinates3(double startCurvature, double endCurvature, double tangentArcStart,
                          double tangentArcEnd, MXClothoid::Coordinates c1, MXClothoid::Coordinates c2,
                          double length, int npts){
        int branchAngleDirection = 0;
        if (startCurvature < endCurvature){
            //第一缓和曲线
            //计算角度差
            double branch_angle = MXClothoid::calcBranchAngle(tangentArcStart, tangentArcEnd);
            int tmp = MXClothoid::convertCurvatureToProtocol(startCurvature);
            branchAngleDirection = ANGLE_DIR(branch_angle, tmp);  //0: right, 1, left
        }else{
            //第二缓和曲线
            //计算角度差
            double branch_angle = MXClothoid::calcBranchAngle(tangentArcEnd, tangentArcStart);
            int tmp = MXClothoid::convertCurvatureToProtocol(endCurvature);
            branchAngleDirection = ANGLE_DIR(branch_angle, tmp);  //0: right, 1, left
        }
        int sign = ANGLE_SIGN(branchAngleDirection);
        //计算地图缩放比例
        double coefficient = getCoefficient(startCurvature, endCurvature, length, sign, c1, c2);

        std::vector<MXClothoid::Coordinates> veccoord;
        std::ofstream outFile;
        outFile.open("../data.csv", std::ios::app|std::ios::out);
        //计算独立坐标相对偏转角度,用于坐标转换
        double diff = MXClothoid::calcClothoidBranchAngle(startCurvature, endCurvature, length, sign);
        double startAngle = (90.0f - ((double)tangentArcStart * 360.0f / 65535)) * PI / 180.0f;
        double endAngle_tmp = (90.0f - ((double)tangentArcEnd * 360.0f / 65535)) * PI / 180.0f;
        double endAngle = startAngle + diff;
        double angle = 0.0;
        if(startCurvature < endCurvature){
            angle = MXClothoid::calcAngleByStartTangentArc(startCurvature, endCurvature, length, sign, startAngle);
        }else{
            angle = MXClothoid::calcAngleByEndTangentArc(startCurvature, endCurvature, length, sign, endAngle);
        }

        double asquare = fabs(length / (endCurvature - startCurvature));
        double La = asquare * startCurvature;
        double Le = asquare * endCurvature;
        double averageLength =  length / npts;
        for (int i = 0; i <= npts; i++) {
            double a = 0.0;
            double b = 0.0;
            if(La > Le){
                a = La;
                b = La - averageLength * i;
            }else{
                a = La;
                b = La + averageLength * i;
            }
            double tempx = MXClothoid::integral(MXClothoid::cosFunction, a, b, asquare);
            double tempy = MXClothoid::integral(MXClothoid::sinFunction, a, b, asquare);

            MXClothoid::Coordinates coord;
            coord.x  = ((tempx * cos(angle) - (sign * tempy) * sin(angle)));
            coord.y  = ((tempx * sin(angle) + (sign * tempy) * cos(angle)));

            veccoord.push_back(coord);
        }
        for (int i = 0; i < veccoord.size(); i++) {
            double x = c1.x + veccoord[i].x * coefficient;
            double y = c1.y + veccoord[i].y * coefficient;
            outFile << std::setprecision(15)<< x << ',' << std::setprecision(15)<< y << std::endl;
            std::cout<<i <<","<< std::setprecision(15) << x <<","<<std::setprecision(15) << y <<",11"<< std::endl;
        }
        outFile.close();
    }

    //验证DB
    void calcCoordinates(int startCurvature, int endCurvature, double length, double tangentArcStart,
            double tangentArcEnd, int npts, MXClothoid::Coordinates c1, MXClothoid::Coordinates c2){
        if(startCurvature == endCurvature){
            std::cout<<"起点曲率不能等于终点曲率!!!"<<std::endl;
            return;
        }
        double sgStartCurvature = fabs(MXClothoid::getcurvature(startCurvature));
        double sgEndCurvature = fabs(MXClothoid::getcurvature(endCurvature));
#if 0
        std::cout<<"-----------------通过起点终点坐标计算----------------------"<<std::endl;
        calcCoordinates1(sgStartCurvature, sgEndCurvature, tangentArcStart, tangentArcEnd,
                         c1, c2, length, npts);
#endif
        std::cout<<"-----------------通过切线方位角度计算----------------------"<<std::endl;
#if 0
        calcCoordinates2(sgStartCurvature, sgEndCurvature, tangentArcStart, tangentArcEnd,
                c1, c2, length, npts);
#endif

#if 1
        calcCoordinates3(sgStartCurvature, sgEndCurvature, tangentArcStart, tangentArcEnd,
                c1, c2, length, npts);
#endif
        std::cout<<"========================================================"<<std::endl;
    }

    void test1(){
        int startX = 1388724896;
        int startY = 475678560;
        int startdx = -1;
        int startdy = -1;

        int endX = 1388730144;
        int endY = 475679744;
        int enddx = -22;
        int enddy = -6;

        int node = 10;                   //分成等份
        int startCurvature = 780;       //开始曲率
        int endCurvature = 215;         //结束曲率
        double tangentArcStart = 13392;
        double tangentArcEnd = 12611;
        double length = 34;             //长度 m
        MXClothoid::Coordinates coord1 = MXClothoid::calcOriginCoordinates(startX, startY, startdx, startdy);
        MXClothoid::Coordinates coord2 = MXClothoid::calcOriginCoordinates(endX, endY, enddx, enddy);
        calcCoordinates(startCurvature, endCurvature, length, tangentArcStart, tangentArcEnd, node, coord1, coord2);
    }

    void test2(){
        int startX = 1388730144;
        int startY = 475679744;
        int startdx = -22;
        int startdy = -6;

        int endX = 1388730336;
        int endY = 475684896;
        int enddx = 0;
        int enddy = -2;

        int node = 10;                   //分成等份
        int startCurvature = 160;       //开始曲率
        int endCurvature = 154;         //结束曲率
        double tangentArcStart = 12611;
        double tangentArcEnd = 55159;
        double length = 60;             //长度 m
        MXClothoid::Coordinates coord1 = MXClothoid::calcOriginCoordinates(startX, startY, startdx, startdy);
        MXClothoid::Coordinates coord2 = MXClothoid::calcOriginCoordinates(endX, endY, enddx, enddy);
        calcCoordinates(startCurvature, endCurvature, length, tangentArcStart, tangentArcEnd, node, coord1, coord2);
    }

    void test3(){
        int startX = 1388730336;
        int startY = 475684896;
        int startdx = 0;
        int startdy = -2;

        int endX = 1388727264;
        int endY = 475685920;
        int enddx = -1;
        int enddy = -1;

        int node = 10;                   //分成等份
        int startCurvature = 178;       //开始曲率
        int endCurvature = 798;         //结束曲率
        double tangentArcStart = 55159;
        double tangentArcEnd = 53346;
        double length = 24;             //长度 m
        MXClothoid::Coordinates coord1 = MXClothoid::calcOriginCoordinates(startX, startY, startdx, startdy);
        MXClothoid::Coordinates coord2 = MXClothoid::calcOriginCoordinates(endX, endY, enddx, enddy);
        calcCoordinates(startCurvature, endCurvature, length, tangentArcStart, tangentArcEnd, node, coord1, coord2);
    }

    void test4(){
        int startX = 1388539680;
        int startY = 475614112;
        int startdx = -1;
        int startdy = -1;

        int endX = 1388539616;
        int endY = 475616448;
        int enddx = 2;
        int enddy = -1;

       // int sign = -1;                   //右转方位角减小 -1 左转方位角减小 -1
        int node = 10;                   //分成等份
        int startCurvature = 262;       //开始曲率
        int endCurvature = 884;         //结束曲率
        double tangentArcStart = 64023;
        double tangentArcEnd = 3682;
        double length = 21;             //长度 m
        MXClothoid::Coordinates coord1 = MXClothoid::calcOriginCoordinates(startX, startY, startdx, startdy);
        MXClothoid::Coordinates coord2 = MXClothoid::calcOriginCoordinates(endX, endY, enddx, enddy);
        calcCoordinates(startCurvature, endCurvature, length, tangentArcStart, tangentArcEnd, node, coord1, coord2);
    }

    void test5(){
        int startX = 1388539616;
        int startY = 475616448;
        int startdx = 2;
        int startdy = -1;

        int endX = 1388541280;
        int endY = 475617920;
        int enddx = -1;
        int enddy = -1;

        //int sign = 1;                   //右转方位角减小 -1 左转方位角减小 -1
        int node = 10;                   //分成等份
        int startCurvature = 892;       //开始曲率
        int endCurvature = 307;         //结束曲率
        double tangentArcStart = 3682;
        double tangentArcEnd = 8700;
        double length = 18;             //长度 m
        MXClothoid::Coordinates coord1 = MXClothoid::calcOriginCoordinates(startX, startY, startdx, startdy);
        MXClothoid::Coordinates coord2 = MXClothoid::calcOriginCoordinates(endX, endY, enddx, enddy);
        calcCoordinates(startCurvature, endCurvature, length, tangentArcStart, tangentArcEnd, node, coord1, coord2);
    }

    void test6(){
        int startX = 1388712320;
        int startY = 475156416;
        int startdx = -1;
        int startdy = -1;

        int endX = 1388713472;
        int endY = 475154656;
        int enddx = -4;
        int enddy = 4;

        int node = 10;                   //分成等份
        int startCurvature = 205;       //开始曲率
        int endCurvature = 118;         //结束曲率
        double tangentArcStart = 31145;
        double tangentArcEnd = 22789;
        double length = 17;             //长度 m
        MXClothoid::Coordinates coord1 = MXClothoid::calcOriginCoordinates(startX, startY, startdx, startdy);
        MXClothoid::Coordinates coord2 = MXClothoid::calcOriginCoordinates(endX, endY, enddx, enddy);
        calcCoordinates(startCurvature, endCurvature, length, tangentArcStart, tangentArcEnd, node, coord1, coord2);
    }

    void test7(){
        int startX = 1388716192;
        int startY = 475155008;
        int startdx = -4;
        int startdy = -1;

        int endX = 1388716704;
        int endY = 475156704;
        int enddx = -1;
        int enddy = -1;

        int node = 10;                   //分成等份
        int startCurvature = 116;       //开始曲率
        int endCurvature = 236;         //结束曲率
        double tangentArcStart = 8191;
        double tangentArcEnd = 65409;
        double length = 17;             //长度 m
        MXClothoid::Coordinates coord1 = MXClothoid::calcOriginCoordinates(startX, startY, startdx, startdy);
        MXClothoid::Coordinates coord2 = MXClothoid::calcOriginCoordinates(endX, endY, enddx, enddy);
        calcCoordinates(startCurvature, endCurvature, length, tangentArcStart, tangentArcEnd, node, coord1, coord2);
    }

    void test8(){
        int startX = 1388712320;
        int startY = 475156416;
        int startdx = -1;
        int startdy = -1;

        int endX = 1388713152;
        int endY = 475157568;
        int enddx = -6;
        int enddy = -7;

        int node = 10;                   //分成等份
        int startCurvature = 696;       //开始曲率
        int endCurvature = 897;         //结束曲率
        double tangentArcStart = 4082;
        double tangentArcEnd = 7445;
        double length = 9;             //长度 m
        MXClothoid::Coordinates coord1 = MXClothoid::calcOriginCoordinates(startX, startY, startdx, startdy);
        MXClothoid::Coordinates coord2 = MXClothoid::calcOriginCoordinates(endX, endY, enddx, enddy);
        calcCoordinates(startCurvature, endCurvature, length, tangentArcStart, tangentArcEnd, node, coord1, coord2);
    }

    void test9(){
        int startX = 1388713152;
        int startY = 475157568;
        int startdx = -6;
        int startdy = -7;

        int endX = 1388716160;
        int endY = 475157408;
        int enddx = -3;
        int enddy = -1;

        int node = 10;                   //分成等份
        int startCurvature = 900;       //开始曲率
        int endCurvature = 897;         //结束曲率
        double tangentArcStart = 7445;
        double tangentArcEnd = 24736;
        double length = 25;             //长度 m
        MXClothoid::Coordinates coord1 = MXClothoid::calcOriginCoordinates(startX, startY, startdx, startdy);
        MXClothoid::Coordinates coord2 = MXClothoid::calcOriginCoordinates(endX, endY, enddx, enddy);
        calcCoordinates(startCurvature, endCurvature, length, tangentArcStart, tangentArcEnd, node, coord1, coord2);
    }
}

#endif //MXCLOTHOID_TEST_H
