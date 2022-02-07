#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>


typedef struct {
	float * pData;
	short int rows;
	short int columns;

}t_matrix;

#define fsize 4
#define pi 3.14159265358979

#define pos3(y,x) x+(y*3)
#define pos4(y,x) x+(y*4)

void errorIdler(){
	printf("error occured");
	while (1);
}

t_matrix createMatrix(int rows,int columns) {
	t_matrix curMat = { malloc(rows * columns * fsize),rows,columns };
	memset(curMat.pData, 0, rows * columns * fsize);
	return(curMat);
}  

void freeMat(t_matrix * mat) {
	free(mat->pData);
}

uint8_t checkSize(t_matrix* pMat, short int desiredRows, short int desiredColumns,char * matName) {

	if ((pMat->rows == desiredRows)&&(pMat->columns == desiredColumns)) {
		return(0);
	}
	else {
		printf("%s matrix is of incorrect size, expected %u x %u,got %u x %u\r\n", matName, desiredRows, desiredColumns, pMat->rows, pMat->columns);
		errorIdler();
		return(1);
	}


}

void setToZero(t_matrix* mat) {
	memset(mat->pData, 0, mat->columns * mat->rows * fsize);
}

void setToIdentity(t_matrix* mat) {
	memset(mat->pData, 0, mat->columns * mat->rows * fsize);
	short int len = mat->columns;
	if (mat->rows < len) { len = mat->rows; }
	for (short int i = 0; i < len; i++) {
		mat->pData[i + (i * mat->columns)] = 1;
	}
}

void printMatrix(t_matrix * pMat,char * name) {
	printf("%s = [\r\n", name);
	for (short int y = 0; y < pMat->rows; y++) {
		for (short int x = 0; x < pMat->columns; x++) {
			printf("%5.4f", pMat->pData[x + (y * pMat->columns)]);
			if (x != pMat->rows-1) { printf(","); }
		}
		printf(";\r\n");
	}
	printf("]\r\n");
}

void addToMatrix(t_matrix* dest, t_matrix* src) {
	int size = dest->columns * dest->rows,srcSize= src->columns * src->rows;
	if (srcSize < size) { size = srcSize; }
	for (int i = 0; i < size; i++) {
		dest->pData[i] += src->pData[i];
	}
}

float dotProduct(float * a,float * b,short int length,short int incOnB){
	//takes the dot product of the two arrays
	//increment on b is the array increment per element summed,
	//this allows for directly reading form a 2d matrix
	float sum=0;
	for (short int i = 0; i < length; i++) {
		sum += a[i] * b[i * incOnB];
		
	}
	return(sum);
}

void multiMat(t_matrix * destMat,t_matrix* matA,t_matrix *matB) {
	//multiplies matrix a by b and returns the result
	if (matA->columns != matB->rows) {
		//printf("Cannot multiply, incorrect size");
		errorIdler();
	}
	short int rows= matA->rows, cols= matB->columns;
	float * result=malloc(rows * cols*fsize);
	
	for (short int y = 0; y <rows; y++) {
		for (short int x = 0; x < cols; x++) {
			result[x+ (y*cols)]= dotProduct(&matA->pData[y*matA->columns],&matB->pData[x],matA->columns,cols);
			//gets the dot product and writes it
		}
	}
	
	memcpy(destMat->pData, result, rows * cols * fsize);
	free(result);
}

void multiMatScalar(t_matrix* matrixOut,t_matrix* matrixIn,float scalar) {
	//multiplies matrix a by a scalar and writes it to the output matrix
	short int rows = matrixIn->rows, cols = matrixIn->columns;
	checkSize(matrixOut, rows, cols, "scalar multiplication result");
	for (short int y = 0; y < rows; y++) {
		for (short int x = 0; x < cols; x++) {
			matrixOut->pData[x + (y*cols)]= matrixIn->pData[x + (y * cols)]*scalar;
		}
	}

}

void multiMatScalarAccumulative(t_matrix* matrixOut, t_matrix* matrixIn, float scalar) {
	//multiplies matrix a by a scalar and adds it to the output matrix
	short int rows = matrixIn->rows, cols = matrixIn->columns;
	checkSize(matrixOut, rows, cols, "scalar multiplication result");
	for (short int y = 0; y < rows; y++) {
		for (short int x = 0; x < cols; x++) {
			matrixOut->pData[x + (y * cols)] += matrixIn->pData[x + (y * cols)] * scalar;
		}
	}

}

float det2x2(t_matrix* matIn, short int avoidRow, short int avoidCol) {
//does a 2x2 subsection of a 3x3 matrix
	checkSize(matIn, 3, 3, "for determinant");
	float sum=0;
	int curIndex = 0;
	for (short int y = 0; y < 3; y++) {
		if (y != avoidRow) {
			for (short int x = 0; x < 3; x++) {
				if (x != avoidCol) {
					sum += matIn->pData[pos3(y, x)] * ((0b11 ^ curIndex) ?  -1: 1);//becomes ind:val= 0:1 1:-1 2:-1 3:1
					curIndex++;
				}
			}
		}
	}
	return(sum);
}

t_matrix adjoint(t_matrix* matIn) {
	checkSize(matIn, 3, 3, "mat for adjoint");
	t_matrix adjMat=createMatrix(3, 3);
	short signed int sign = 1;
	for (short int y = 0; y < matIn->rows; y++) {
		for (short int x = 0; x < matIn->columns; x++) {
			adjMat.pData[pos3(y,x)]=det2x2(matIn, y, x)*sign;
			sign *= -1;
			//add the determinant of all elements in other rows and flip the sign
		}
	}
	return(adjMat);
}


t_matrix dupeMat(t_matrix* matIn) {
	t_matrix matOut = createMatrix(matIn->rows, matIn->columns);
	memcpy(matOut.pData, matIn->pData, matIn->columns * matIn->rows * fsize);
	return(matOut);
}

void transpose(t_matrix* matOut, t_matrix* matIn) {
	checkSize(matOut, matIn->columns, matIn->rows, "transpose");
	for (short int y = 0; y < matIn->rows; y++) {
		for (short int x = 0; x < matIn->columns; x++) {
			matOut->pData[y + (matIn->columns * x)] = matIn->pData[x + (matIn->columns * y)];
		}
	}
}

void insertRotMat(t_matrix* tfMat, t_matrix* rotMat) {
	//Inserts a rotation matrix into the upper left corner of a transform matrix
	checkSize(tfMat, 4, 4, "transform");
	checkSize(rotMat, 3, 3, "rotation");

	memcpy(tfMat->pData, rotMat->pData, fsize * 3);
	memcpy(tfMat->pData + (fsize * pos4(1,0)), rotMat->pData + (fsize * pos3(1,0)), fsize * 3);
	memcpy(tfMat->pData + (fsize * pos4(2,0)), rotMat->pData + (fsize * pos3(2,0)), fsize * 3);
}

void extractRotMat(t_matrix* rotMat, t_matrix* tfMat) {
	//Inserts a rotation matrix into the upper left corner of a transform matrix
	checkSize(tfMat, 4, 4, "transform");
	checkSize(rotMat, 3, 3, "rotation");

	memcpy(rotMat->pData, tfMat->pData, fsize * 3);
	memcpy(rotMat->pData + (fsize * pos3(1, 0)), tfMat->pData + (fsize * pos4(1, 0)), fsize * 3);
	memcpy(rotMat->pData + (fsize * pos3(2, 0)), tfMat->pData + (fsize * pos4(2, 0)), fsize * 3);
}

void insertXYZtoTF(t_matrix* tfMat, float* xyz) {
	checkSize(tfMat, 4, 4, "transform");
	tfMat->pData[pos4(0, 3)] = xyz[0];
	tfMat->pData[pos4(1, 3)] = xyz[1];
	tfMat->pData[pos4(2, 3)] = xyz[2];
}

void extractXYZfromTF(float* xyz, t_matrix* tfMat) {
	checkSize(tfMat, 4, 4, "transform");
	xyz[0]=tfMat->pData[pos4(0, 3)];
	xyz[1]=tfMat->pData[pos4(1, 3)];
	xyz[2]=tfMat->pData[pos4(2, 3)];

}

void vecToSO3(t_matrix * skewOut,float* pVec) {
	//takes angular velocity 3 element vector and returns skew matrix
	
	checkSize(skewOut, 3, 3, "skew");
	float matValues[9] = {
		0,			-pVec[2],	pVec[1],
		pVec[2],	0,			-pVec[0],
		-pVec[1],	pVec[0],	0};
	memcpy(skewOut->pData,matValues,9*fsize);
	
	
}


void so3ToVec(float* vecOut, t_matrix* skewIn) {
	checkSize(skewIn, 3, 3, "skew");
	vecOut[0] = skewIn->pData[pos3(2, 1)];
	vecOut[1] = skewIn->pData[pos3(0, 2)];
	vecOut[2] = skewIn->pData[pos3(1, 0)];
}


void vecToSE3(t_matrix* se3Out,float * pVec) {
	//spatial velocity 6 element vector to 4x4 se3 matrix
	checkSize(se3Out, 4, 4, "se3");
		float matValues[16] = {
			0,0,0,pVec[3],
			0,0,0,pVec[4],
			0,0,0,pVec[5],
			0,0,0,0};
		//copy in zeros and linear velocity
		memcpy(se3Out->pData, matValues, 16 * fsize);
		t_matrix skewMat = createMatrix(3, 3);
		vecToSO3(&skewMat, pVec);
		insertRotMat(se3Out, &skewMat);
}

float normalise(float * pVec,short int length) {
	//gets absolute magnitude of a 2d or 3d vector
	float curSum=0;
	for (short int i = 0; i < length; i++) {
		curSum += (pVec[i] * pVec[i]);
	}
	return(sqrtf(curSum));
}

float axisAng3(float * rotAxis,float* pVec) {
	//takes rotation 3 element vector and returns unit rot axis and rotation angle
	float theta = normalise(pVec,3);
	if (rotAxis != NULL) {
		rotAxis[0] = pVec[0] / theta;
		rotAxis[1] = pVec[1] / theta;
		rotAxis[2] = pVec[2] / theta;
	}
	return(theta);
}

#define threshold 10^-7
void roundNearZeroToZero(t_matrix* pMat) {
	for (short int y = 0; y < pMat->rows; y++) {
		for (short int x = 0; x < pMat->columns; x++) {
			if ((abs(pMat->pData)) < threshold) { pMat->pData = 0; }
		}
	}

}


void rotMatrixExp(t_matrix* rotOut, t_matrix* so3RotIn) {
	checkSize(rotOut, 3, 3, "Exponential rotation");
	checkSize(rotOut, 3, 3, "so3 rotation");
	
		float angleVector[3];
		so3ToVec(angleVector, so3RotIn);
		float angleAboutUnit = axisAng3(NULL, angleVector);
		t_matrix unitRotMat = createMatrix(3, 3);
		multiMatScalar(&unitRotMat, so3RotIn, 1 / angleAboutUnit);
		//divide the rotation matrix by theta to get the unit rotation matrix
		setToIdentity(rotOut);
		//sum identity matrix,
		//sin(angle) * unitRot (1st order terms)
		multiMatScalarAccumulative(rotOut, &unitRotMat, sinf(angleAboutUnit));

		// and (1-cos(angle))*unitRot*unitRot (2nd order terms)

		multiMat(&unitRotMat, &unitRotMat, &unitRotMat);
		multiMatScalarAccumulative(rotOut, &unitRotMat, 1 - cosf(angleAboutUnit));
	
}

void tfMatrixExp(t_matrix* Tse3out, t_matrix* se3PosIn) {
	checkSize(Tse3out, 4, 4, "Screw axis from exponential");
	checkSize(se3PosIn, 4, 4, "exponential transform");
	float posToTF[3];
	t_matrix posRot = createMatrix(3, 3);
	extractRotMat(&posRot, se3PosIn);
	//printMatrix(&posRot, "Rotation matrix");
	if (normalise(posRot.pData, 9) < 0.000001) {
		setToIdentity(Tse3out);//avoids divide by zero errors by just setting it to a identity matrix
		Tse3out->pData[pos4(0, 3)] = se3PosIn->pData[pos4(0, 3)];
		Tse3out->pData[pos4(1, 3)] = se3PosIn->pData[pos4(1, 3)];
		Tse3out->pData[pos4(2, 3)] = se3PosIn->pData[pos4(2, 3)];
		return;
	}
	else {
		//extract the rotation matrix from the transformation matrix coming in
		so3ToVec(posToTF, &posRot);
		//printf("posToTf is %f, %f, %f \r\n", posToTF[0], posToTF[1], posToTF[2]);
		//get the magnitude of rotation about those three axis
		float curRotAngle = axisAng3(NULL, posToTF);
		//printf("cur rot angle (theta) is %E \r\n", curRotAngle);

		t_matrix unitRotAxis = createMatrix(3, 3);

		multiMatScalar(&unitRotAxis, &posRot, 1 / curRotAngle);
		//printMatrix(&unitRotAxis, "Unit rotation matrix");
		//divide the rotation matrix by that angle to get a unit rotation matrix

			//convert rotation to exponential form and insert
		t_matrix expoRot = createMatrix(3, 3);
		rotMatrixExp(&expoRot, &posRot);
		insertRotMat(Tse3out, &expoRot);
		//printMatrix(&expoRot, "exponential rotation matrix");
		t_matrix accumulatorMat = createMatrix(3, 3);

		//create identity matrix times the angle of rotation and put it into a tempoary matrix
		accumulatorMat.pData[0] = curRotAngle;
		accumulatorMat.pData[pos3(1, 1)] = curRotAngle;
		accumulatorMat.pData[pos3(2, 2)] = curRotAngle;

		multiMatScalarAccumulative(&accumulatorMat, &unitRotAxis, 1.0f - cosf(curRotAngle));
		//first order terms are added (unit rotation axis matrix times 1-cos(angle)


		t_matrix unitRotScaled = createMatrix(3, 3);

		multiMatScalar(&unitRotScaled, &unitRotAxis, curRotAngle - sinf(curRotAngle));
		multiMat(&unitRotScaled, &unitRotScaled, &unitRotAxis);
		//second order terms (rotation axis squared times angle-sin(angle)) are added
		addToMatrix(&accumulatorMat, &unitRotScaled);
		t_matrix xyzVec = createMatrix(3, 1);

		extractXYZfromTF(xyzVec.pData, se3PosIn);

		multiMat(&xyzVec, &accumulatorMat, &xyzVec);
		multiMatScalar(&xyzVec, &xyzVec, (1.0f / curRotAngle));

		//printMatrix(&xyzVec, "xyz vector");
		insertXYZtoTF(Tse3out, xyzVec.pData);
		
		Tse3out->pData[15] = 1;
		freeMat(&expoRot);
		freeMat(&accumulatorMat);
		freeMat(&unitRotAxis);
		//printMatrix(Tse3out, "TseOut");
		//garbage collection
	}
	freeMat(&posRot);
	
	
		
}

void FKinSpace(t_matrix * toolPos,t_matrix* pHomeTF, t_matrix * Slist, float* thetaList,short int numThetas) {
	t_matrix expT=createMatrix(4,4);
	memcpy(expT.pData,pHomeTF->pData,(pHomeTF->columns* pHomeTF->rows)*fsize);
	//copy the initial home configuration into the working matrix
	t_matrix screwAxisInExp=createMatrix(4,4);
	t_matrix SscaledToTheta = createMatrix(Slist[0].rows, Slist[0].columns);
	t_matrix expJointTf = createMatrix(4, 4);
	for(signed int i = numThetas-1; i >= 0; i--) {
		//printf("\r\n i=%u\r\n", i);
		//printMatrix(&expT, "T");

		//working from last joint (end effector) to first (nearest base)
		multiMatScalar(&SscaledToTheta,&Slist[i],thetaList[i]);
		//printf("theta is %f\r\n", thetaList[i]);
		//printMatrix(&Slist[i], "Screw axis");
		//printMatrix(&SscaledToTheta, "Screw axis scaled");

		vecToSE3(&screwAxisInExp, SscaledToTheta.pData);
		//printMatrix(&screwAxisInExp, "Screw axis in exp");
		tfMatrixExp(&expJointTf, &screwAxisInExp);
		//printMatrix(&expJointTf, "joint tf exponential");
		multiMat( & expT, & expJointTf, & expT);


		//the transform is equal to the exponential of the joint screw axis times the joint coordinates multiplied by the previous transform

	}
	
	memcpy(toolPos->pData,expT.pData,16*fsize);
	freeMat(&expT);
	freeMat(&screwAxisInExp);
	freeMat(&SscaledToTheta);
	freeMat(&expJointTf);
	
}

void fastTFinverse(t_matrix * matOut, t_matrix* matIn) {

	//utilises tf format to take the inverse faster
	checkSize(matOut, 4, 4, "transfer function inverse");
	checkSize(matIn, 4, 4, "transfer function to be inverted");
	memset(matOut->pData, 0, fsize * 16);
	t_matrix rotMat = createMatrix(3, 3);
	extractRotMat(&rotMat, matIn);
	
	t_matrix rotMatTranspose = createMatrix(3, 3);
	transpose(&rotMatTranspose,&rotMat);
	freeMat(&rotMat);
	insertRotMat(matOut,&rotMatTranspose);
	multiMatScalar(&rotMatTranspose, &rotMatTranspose, -1);
	t_matrix xyzVec=createMatrix(3,1);
	extractXYZfromTF(xyzVec.pData, matIn);
	multiMat(&xyzVec,&rotMatTranspose,&xyzVec);
	insertXYZtoTF(matOut, &xyzVec);
	matOut->pData[15] = 1;
	freeMat(&rotMatTranspose);

}

t_matrix IKinSpace(t_matrix* desiredTF,t_matrix * Slist, t_matrix* homeTF,t_matrix * initAngles,short int numAngles,float linError,float rotError) {
	t_matrix prevJointAngles = dupeMat(initAngles);
	t_matrix vS;
	short int iterations = 0;
	do {


		iterations++;
	} while ((normalise(vS.pData,3))>rotError || normalise(vS.pData+(fsize*3),3)>linError);

}

int main(char* args) {

	float testHomeTFdata[] = {
	-1, 0, 0, 0 ,
	0, 1, 0, 6,
	0, 0, -1, 2,
	0, 0, 0, 1
	};
	t_matrix testHomeTF =
	{ testHomeTFdata,4,4 };
	printMatrix(&testHomeTF, "home transform mat");
	float SlistData[3][6] = {
		{0,0,1,4,0,0},
		{0,0,0,0,1,0},
		{0,0,-1,-6,0,-0.1}
	};

	t_matrix Slist[6] = {
		{SlistData[0],6,1},
		{SlistData[1],6,1},
		{SlistData[2],6,1}
	};
	float thetaList [3] = {pi / 2, 3, pi};
	
	t_matrix testOutput = createMatrix(4, 4);
	FKinSpace(&testOutput, &testHomeTF, Slist, thetaList, 3);
	
	//printMatrix(&testMatB, "B");
	printMatrix(&testOutput, "Result");
	
	return(0);
}