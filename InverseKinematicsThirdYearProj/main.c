#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

typedef struct {
	float * pData;
	short int rows;
	short int columns;

}t_matrix;

#define fsize 4
#define pi 3.14159265358979

#define pos3(y,x) x+(y*3)
#define pos4(y,x) x+(y*4)
//sets up 2d indexing for 3x3 and 4x4 matrices
#define jacoR 6
#define jacoC 6
//defines size of jacobian matrix


void errorIdler(){
	printf("error occured");
	while (1);
}

t_matrix createMatrix(int rows,int columns) {
	t_matrix curMat = { malloc(rows * columns * fsize),rows,columns };
	memset(curMat.pData, 0, rows * columns * fsize);
	return(curMat);
}  

#define freeMat(m) free(m.pData)

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
	float * result = malloc(rows * cols * sizeof(float));
	
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
	
	for (short int y = 0; y < 3; y++) {
		if (y != avoidRow) {
			for (short int x = 0; x < 3; x++) {
				if (x != avoidCol) {
					sum += matIn->pData[pos3(y, x)] * ((y+x)%2)?-1:1;//becomes ind:val= 0:1 1:-1 2:-1 3:1
					
				}
			}
		}
	}
	return(sum);
}

void extractSubMatrix(t_matrix* matOut, t_matrix* matIn,short int rowToAvoid,short int colToAvoid) {
	checkSize(matOut, matIn->rows - 1, matIn->columns - 1, "sub ");
	short int curIndex = 0;
	for (short int y = 0; y < matIn->rows; y++) {
		if (y != rowToAvoid) {
			for (short int x = 0; x < matIn->columns; x++) {
				if (x != colToAvoid) {
					matOut->pData[curIndex] = matIn->pData[x + (matIn->columns * y)];
				}

			}
		}
	}
}

float detNxN(t_matrix* matIn) {
	float sum = 0;
	
	short int rows = matIn->rows,cols=matIn->columns;
	//recurses down to a 2x2 matrix
	t_matrix subMat = createMatrix(rows - 1, cols - 1);
	for (short int y = 0; y < rows; y++) {
		for (short int x = 0; x < cols; x++) {
			if (rows == 2 && cols == 2) {
				sum = matIn->pData[0] + matIn->pData[3] - (matIn->pData[1]+matIn->pData[2]);
			}
			else {
				extractSubMatrix(&subMat, matIn, y, x);
				sum = matIn->pData[x + (y * cols)] * (((y + x) % 2) ? 1 : -1) * detNxN(&subMat);
			}
		}
	}	

	freeMat(subMat);
	return(sum);

}



void cofactor(t_matrix * matOut,t_matrix* matIn) {
	checkSize(matOut, matIn->rows, matIn->columns, "mat for adjoint");
	float * matBuf = malloc(fsize * matIn->rows * matIn->columns);
	short signed int sign = 1;
	t_matrix subMat = createMatrix(matIn->rows - 1, matIn->columns - 1);
	for (short int y = 0; y < matIn->rows; y++) {
		for (short int x = 0; x < matIn->columns; x++) {
			extractSubMatrix(&subMat, matIn, y, x);
			matBuf[(matIn->columns*y) +x]=detNxN(&subMat)* (((y + x) % 2)) ? -1 : 1;
			
			//add the determinant of all elements in other rows
		}
	}

	memcpy(matOut->pData, matBuf, fsize * matIn->rows * matIn->columns);
	free(matBuf);
	free(&subMat);
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
/*
* old adjoint
void adjoint(t_matrix* matOut, t_matrix* matIn) {
	checkSize(matOut, matIn->columns, matIn->rows, "adjoint");
	t_matrix cofactorMat = createMatrix(matIn->rows, matIn->columns);
	cofactor(&cofactorMat, matIn);
	transpose(matOut, &cofactorMat);
	freeMat(cofactorMat);
}
*/

/*copies one size matrix variables to another whilst keeping rows intact
align x and y trim the position of the source matrix within the destination that is written to
*/
void copyDiffSizeMatToMat(t_matrix* dest, t_matrix* src,short int alignY,short int alignX) {
	//assumes destination is always large enough to fit source
	short int curRow = 0,curCol=0;
	for (short int y = 0; y < src->rows; y++) {
		memcpy(&(dest->pData[(alignX + (dest->columns * (alignY + y)))]), &(src->pData[y* src->columns]), fsize * src->columns);
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
	vecOut[0] = skewIn->pData[(2*skewIn->columns)+1];
	vecOut[1] = skewIn->pData[ 2];
	vecOut[2] = skewIn->pData[skewIn->columns +0];
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
		freeMat(skewMat);
}

void se3ToVec(float* pVec, t_matrix* se3In) {
	//printMatrix(se3In, "se3in");
	//gets 6 element velocity vector from se3
	checkSize(se3In, 4, 4, "se3 mat to have vector found");
	pVec[0] = se3In->pData[(2 * se3In->columns)+1];
	pVec[1] = se3In->pData[2];
	pVec[2] = se3In->pData[se3In->columns];
	pVec[3] = se3In->pData[3];//x
	pVec[4] = se3In->pData[se3In->columns +3];//y
	pVec[5] = se3In->pData[(se3In->columns*2)+3];//z

}

float normalise(float * pVec,short int length) {
	//gets absolute magnitude of a 2d or 3d vector
	float curSum=0;
	for (short int i = 0; i < length; i++) {
		curSum += (pVec[i] * pVec[i]);
	}
	return(sqrtf(curSum));
}

float normDiffHomogenousTFrot(t_matrix* matA, t_matrix* matB) {
	float sum = 0;
	for (short int y = 0; y <3; y++)
	{
		for (short int x = 0; x <3; x++)
		{
			sum+=pow(matA->pData[pos4(y, x)] - matB->pData[pos4(y, x)], 2);
		}
	}
	return(sqrt(sum));
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
			if (((abs(pMat->pData))) < threshold) { pMat->pData = 0; }
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
	//takes the log of a exponential
	checkSize(Tse3out, 4, 4, "Screw axis");
	checkSize(se3PosIn, 4, 4, "exponential transform input");
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
		freeMat(expoRot);
		freeMat(accumulatorMat);
		freeMat(unitRotAxis);
		//printMatrix(Tse3out, "TseOut");
		//garbage collection
	}
	freeMat(posRot);
	
	
		
}


void adjoint(t_matrix* matOut, t_matrix* matIn) {
	//puts a se3 transformation into a 6x6 adjoint representation
	checkSize(matIn, 4, 4, "to have adjoint found");
	checkSize(matOut, 6, 6, "adjoint output");
	//printMatrix(matIn, "mat to be made adjoint");
	memset(matOut->pData, 0, 6 * 6 * fsize);
	t_matrix rotMat = createMatrix(3, 3);
	extractRotMat(&rotMat, matIn);
	copyDiffSizeMatToMat(matOut, &rotMat, 0, 0);
	//copy rotation into upper left quadrant
	copyDiffSizeMatToMat(matOut, &rotMat, 3, 3);
	//copy rotation into bottom right quadrant
	t_matrix xyzVec = createMatrix(3, 1);
	extractXYZfromTF(xyzVec.pData, matIn);
	t_matrix skew = createMatrix(3, 3);
	vecToSO3(&skew, xyzVec.pData);
	multiMat(&skew, &skew, &rotMat);
	//convert the position to a so3 representation then post multiply that by the rotation matrix 
	copyDiffSizeMatToMat(matOut, &skew, 3, 0);
	//and put it into the bottom left quadrant
	freeMat(skew);
	freeMat(xyzVec);
	freeMat(rotMat);
}


float traceDiag(t_matrix * mat) {
	//sum of diagonal elements
	short int len = mat->columns;
	float sum = 0;
	if (mat->rows < len) { len = mat->rows; }
	for (short int i = 0; i < len; i++) {
		sum+=mat->pData[i + (i * mat->columns)];
	}
	return(sum);
}



void rotMatrixLog(t_matrix* so3RotOut, t_matrix* rotIn) {
	//turns a matrix into its exponential
	checkSize(rotIn, 3, 3, "rotation for so3 conv");
	checkSize(so3RotOut, 3, 3, "so3 rotation output");
	float rotAngleCos = (traceDiag(rotIn) - 1) / 2;
	if (rotAngleCos >= 1) {
		//if the exponential angle of the input is more than one the output is zero
		memset(so3RotOut->pData, 0, fsize * 3 * 3);
	}
	else if(rotAngleCos <=-1) {
		
	}
	else {
		float rotAngle = acos(rotAngleCos);
		t_matrix product=createMatrix(3,3);
		transpose(&product, rotIn);
		multiMatScalar(&product, &product, -1);
		addToMatrix(&product, rotIn);
		multiMatScalar(so3RotOut,&product, rotAngle*(1/(2*sinf(rotAngle))));
		freeMat(product);
	}
}

void tfMatrixLog(t_matrix* se3PosOut, t_matrix* Tse3in) {
	//takes a matrix and puts it into exponential coordinates
	checkSize(Tse3in, 4, 4, "Expo tf for log");
	checkSize(se3PosOut, 4, 4, "output of log tf");
	memset(se3PosOut->pData + (fsize * pos4(3, 0)), 0, 4 * fsize);
	t_matrix posRotLog = createMatrix(3, 3);
	t_matrix xyzVals = createMatrix(3, 1);
	extractXYZfromTF(xyzVals.pData, Tse3in);
	extractRotMat(&posRotLog, Tse3in);
	t_matrix posRot = createMatrix(3, 3);
	rotMatrixLog(&posRot, &posRotLog);

	if ((normalise(posRot.pData, 9) < 0.000001)) {//avoids divide by zeroes if log of the rotation of the input is negligable
		memset(se3PosOut->pData, 0, 9 * fsize);
		insertXYZtoTF(se3PosOut, xyzVals.pData);
		return;
	}
	
	//printMatrix(&posRotLog, "rot mat for angle");

	float rotAngle = acosf((traceDiag(&posRotLog) - 1) / 2);

	insertRotMat(se3PosOut, &posRot);
	t_matrix accumlator = createMatrix(3, 3);
	setToIdentity(&accumlator);
	multiMatScalarAccumulative(&accumlator, &posRot, -0.5);
	float scalarMult = (1 / rotAngle - (1 / tanf(rotAngle / 2)) / 2)/rotAngle;
	
	multiMat(&posRot, &posRot, &posRot);//square the positional rotation
	multiMatScalarAccumulative(&accumlator, &posRot, scalarMult);
	multiMat(&xyzVals, &accumlator, &xyzVals);//takes the log of the currently exponential coordinates
	insertXYZtoTF(se3PosOut, xyzVals.pData);
	

	
	
}

void FKinSpace(t_matrix * toolPos,t_matrix* pHomeTF, t_matrix * Slist, t_matrix * thetaList,short int numThetas) {
	checkSize(thetaList, numThetas, 1, "list of thetas");
	//printMatrix(thetaList, "forward kin thetas");
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
		multiMatScalar(&SscaledToTheta,&Slist[i],thetaList->pData[i]);
		//printf("theta is %f\r\n", thetaList->pData[i]);
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
	freeMat(expT);
	freeMat(screwAxisInExp);
	freeMat(SscaledToTheta);
	freeMat(expJointTf);
	
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
	freeMat(rotMat);
	insertRotMat(matOut,&rotMatTranspose);
	multiMatScalar(&rotMatTranspose, &rotMatTranspose, -1);
	t_matrix xyzVec=createMatrix(3,1);
	extractXYZfromTF(xyzVec.pData, matIn);
	multiMat(&xyzVec,&rotMatTranspose,&xyzVec);
	
	insertXYZtoTF(matOut, xyzVec.pData);
	matOut->pData[15] = 1;
	freeMat(rotMatTranspose);

}

void jacobianSpace(t_matrix* jacOut, t_matrix* Slist, float* angleLst) {
	//jacobian in space frame - the velocity of the origin
	checkSize(jacOut,6,6,"jacobian");
	t_matrix jacTransposed = createMatrix(6, 6);
	t_matrix jacSingleRow = createMatrix(1, 6);
	//makes addresssing way easier
	for(short int i=0;i<6;i++){
		memcpy(&jacTransposed.pData[6 * i], Slist[i].pData, fsize * 6);
	}
	printMatrix(&jacTransposed, "jac transposed");
	t_matrix Tmat = createMatrix(4, 4);
	t_matrix workingMatrix = createMatrix(4, 4);
	t_matrix adjTmat = createMatrix(6, 6);
	t_matrix workingTwistMat = createMatrix(6, 1);
	t_matrix jointPos = createMatrix(4, 4);
	setToIdentity(&Tmat);
	for (short int i = 0; i < 5; i++) {
		printf("jacobian I=%u\r\n", i);
		printMatrix(&(Slist[i]), "Slist vals");
		multiMatScalar(&workingTwistMat, &(Slist[i]), (angleLst[i]));
		printMatrix(&workingTwistMat, "current twist matrix");
		vecToSE3(&workingMatrix, workingTwistMat.pData);
		printMatrix(&workingMatrix, "se3 rep of sList[i]*angleList[i]");
		tfMatrixExp(&jointPos,&workingMatrix);
		printMatrix(&jointPos, "expo output joint position in jacobian");
		multiMat(&Tmat, &Tmat, &jointPos);
		printMatrix(&Tmat, "T");
		adjoint(&adjTmat, &Tmat);
		printMatrix(&adjTmat, "adjoint of T");
		multiMat(&jacSingleRow, &adjTmat, &Slist[i+1]);
		printMatrix(&jacSingleRow, "adj of T times Slist[i+1]");
		memcpy(&(jacTransposed.pData[6*(i+1)]), jacSingleRow.pData, 6 * fsize);
		printMatrix(&jacTransposed,"jacobian transposed");
		
	}
	printf("finished finding jacobian");
	printMatrix(&Tmat, "T");
	
	transpose(jacOut, &jacTransposed);
	freeMat(jacTransposed);
	freeMat(Tmat);
	freeMat(workingMatrix);
	freeMat(workingTwistMat);
	freeMat(adjTmat);
	freeMat(jacSingleRow);
	freeMat(jointPos);
}

void swapRows(t_matrix* mat,short int row1,short int row2) {
	float buffer;
	short int cols = mat->columns;
	for (short int i = 0; i < cols; i++)
	{
		buffer = mat->pData[i + (cols * row1)];
		mat->pData[i + (cols * row1)]= mat->pData[i + (cols * row2)];
		mat->pData[i + (cols * row2)]=buffer;
	}

}


void swapCols(t_matrix* mat, short int col1, short int col2) {
	float buffer;
	short int rows = mat->rows,cols=mat->columns;
	for (short int i = 0; i < rows; i++)
	{
		buffer = mat->pData[col1 + (cols * i)];
		mat->pData[col1 + (cols * i)] = mat->pData[col2 + (cols * i)];
		mat->pData[col2 + (cols * i)] = buffer;
	}

}

#define posN(Y,X,C) X+(Y*C)
#define inclusionThresh 1E-15
int gsl_print_matrix(const gsl_matrix* m)
{
	int status, n = 0;

	for (size_t i = 0; i < m->size1; i++) {
		for (size_t j = 0; j < m->size2; j++) {
			if ((status = printf( "%g ", gsl_matrix_get(m, i, j))) < 0)
				return -1;
			n += status;
		}

		if ((status = printf( "\n")) < 0)
			return -1;
		n += status;
	}

	return n;
}

void psudoInverseJaco(t_matrix* matOut, t_matrix* matIn) {
	printMatrix(matIn, "input into psudoInverse");
	uint8_t hasBeenTransposedFlag = 0;
	short int rows = matIn->rows, cols = matIn->columns;

	gsl_matrix* matForSVD = gsl_matrix_alloc(rows, cols);

	if (matIn->columns > matIn->rows) {
		//needs to be transposed if its wider than it is tall
		hasBeenTransposedFlag = 1;

		rows = matIn->columns;
		cols = matIn->rows;
		for (uint8_t y = 0; y < rows; y++) {
			for (uint8_t x = 0; x < cols; x++)
			{
				gsl_matrix_set(matForSVD,y,x, matIn->pData[posN(x,y,rows)]);
			}
			
		}
	}
	else {
		//not transposed
		for (uint8_t y = 0; y < rows; y++) {
			for (uint8_t x = 0; x < cols; x++)
			{
				gsl_matrix_set(matForSVD, y, x, matIn->pData[posN(y, x, cols)]);
			}

		}

	}
	//printf("mat for svd = \n");
	//gsl_print_matrix(matForSVD);
	
	gsl_matrix* matV = gsl_matrix_alloc(cols, cols);
	gsl_vector* sVec = gsl_vector_alloc(cols);
	gsl_vector* work = gsl_vector_alloc(cols);
	gsl_linalg_SV_decomp(matForSVD, matV, sVec, work);
	//factorises A into U*S*transpose(V)
	//matForSVD starts as A but is U on output

	//printf("U (matForSVD) = \n");
	//gsl_print_matrix(matForSVD);
	//printf("V = \n");
	//gsl_print_matrix(matV);

	
	gsl_matrix * Sigma_pinv = gsl_matrix_alloc(rows,cols);
	gsl_matrix_set_zero(Sigma_pinv);

	float thresh = inclusionThresh * gsl_vector_max(sVec);

	float tempF;

	for (short int x = 0; x < cols; x++) {
		if (gsl_vector_get(sVec, x) > thresh) {
			tempF = 1.0f / gsl_vector_get(sVec, x);
		}
		else { tempF = 0.0f; }

		gsl_matrix_set(Sigma_pinv, x, x, tempF);
	}

	gsl_matrix * svdPadded = gsl_matrix_alloc(rows, rows);
	gsl_matrix_set_zero(svdPadded);

	for (short int y = 0; y < rows; y++)
	{
		for (short int x = 0; x < cols; x++)
		{
			gsl_matrix_set(svdPadded, y, x, gsl_matrix_get(matForSVD, y, x));
		}
	}
	//printf("svdPadded = \n");
	//gsl_print_matrix(svdPadded);

	gsl_matrix* gslTempMat = gsl_matrix_alloc(cols, rows);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0f, matV, Sigma_pinv, 0.0f, gslTempMat);

	gsl_matrix* A_pinv;
	if (hasBeenTransposedFlag) {
		A_pinv = gsl_matrix_alloc(rows, cols);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0f, svdPadded, gslTempMat, 0.0f, A_pinv);
	}
	else {
		A_pinv = gsl_matrix_alloc(cols, rows);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0f, gslTempMat, svdPadded, 0.0f, A_pinv);
	}
	//printf("A_pinv = \n");
	//gsl_print_matrix(A_pinv);

	for (short int y = 0; y < rows; y++)
	{
		for (short int x = 0; x < cols; x++)
		{
			matOut->pData[posN(y, x, cols)] = gsl_matrix_get(A_pinv, y, x);
		}
	}
	gsl_matrix_free(matForSVD);
	gsl_matrix_free(matV);
	gsl_vector_free(sVec);
	gsl_vector_free(work);
	gsl_matrix_free(Sigma_pinv);
	gsl_matrix_free(svdPadded);
	gsl_matrix_free(gslTempMat);
	gsl_matrix_free(A_pinv);
	printMatrix(matOut, "output of psudoInverse");
	
}

void IKinSpace(float * neededJointAngles,t_matrix* desiredTF,t_matrix * Slist, t_matrix* homeTF,t_matrix * initAngles,float linError,float rotError) {
	//only works for 6 jointed robots
	//printMatrix(initAngles, "recieved initial angles");
	t_matrix angleMat = createMatrix(6, 1);
	memcpy(angleMat.pData, initAngles->pData, fsize * 6);
	t_matrix tsb = createMatrix(4, 4);
	//printMatrix(initAngles, "copied initial angles");
	/*
	* Vs=Adjoint(forwardKin(homeTF)) x se3of(matLog(Inv(forwardKin(homeTF))*desiredPos)
	*/
	t_matrix vS=createMatrix(6,1);
	t_matrix adjtsb = createMatrix(6, 6);
	t_matrix workingMat = createMatrix(4, 4);
	t_matrix curJaco = createMatrix(6, 6);
	t_matrix invOfJaco = createMatrix(6, 6);
	t_matrix matLogOfInvTF = createMatrix(4, 4);
	FKinSpace(&tsb, homeTF, Slist, &angleMat, 6);
	//printMatrix(&tsb, "working forward transform");
	//get current forward kinematic space frame transform
	adjoint(&adjtsb, &tsb);
	//printMatrix(&adjtsb, "adjoint of cur forward kin tf");

	fastTFinverse(&workingMat, &tsb);
	//printMatrix(&workingMat, "IK tf inverse");
	multiMat(&workingMat,&workingMat, desiredTF);
	tfMatrixLog(&matLogOfInvTF, &workingMat);
	//printMatrix(&matLogOfInvTF, "log of forward tf inverse");
	
	se3ToVec(vS.pData, &matLogOfInvTF);
	//printMatrix(&vS, "Vector of prev");
	multiMat(&vS, &adjtsb, &vS);
	//printMatrix(&vS, "vS");
	short int iterations = 0;
	t_matrix angleDeltaMat = createMatrix(6, 1);
#define maxIterations 20
	int* permutations = malloc(sizeof(int) * 7);
	//newton raphsen method
	//conveniently as ive chosen to model a robot with 6 joints our jacobian is square and a psudoinverse isnt needed
	while((((normalise(vS.pData, 3)) > rotError) || (normalise(&(vS.pData[3]), 3) > linError)) && (iterations<maxIterations)) {
		jacobianSpace(&curJaco, Slist, angleMat.pData);
		printMatrix(&curJaco, "current jacobian");
		psudoInverseJaco(&invOfJaco, &curJaco);
		multiMat(&angleDeltaMat,&invOfJaco, &vS);
		printMatrix(&angleDeltaMat, "change in joint angles");
		addToMatrix(&angleMat, &angleDeltaMat);

		FKinSpace(&tsb, homeTF, Slist, &angleMat, 6);
		printMatrix(&tsb, "working forward transform");
		//get current forward kinematic space frame transform
		adjoint(&adjtsb, &tsb);
		//printMatrix(&adjtsb, "adjoint of cur forward kin tf");

		fastTFinverse(&workingMat, &tsb);
		//printMatrix(&workingMat, "IK tf inverse");
		multiMat(&workingMat, &workingMat, desiredTF);
		tfMatrixLog(&matLogOfInvTF, &workingMat);
		//printMatrix(&matLogOfInvTF, "log of forward tf inverse");

		se3ToVec(vS.pData, &matLogOfInvTF);
		//printMatrix(&vS, "Vector of prev");
		multiMat(&vS, &adjtsb, &vS);
		printMatrix(&vS, "vS");
		iterations++; 
	} 

	printMatrix(&tsb,"final forward kinematic matrix");

	if (iterations == maxIterations) { 
		printf("solution could not be found within %u iterations",maxIterations); }
	else { 	
		printf("solution was found at iteration %u\r\n", iterations); 
		multiMatScalar(&tsb, &tsb, -1.0);
		addToMatrix(&tsb, desiredTF);
		printMatrix(&tsb,"forward transform error matrix");
	}

	memcpy(neededJointAngles, angleMat.pData, fsize * 6);
	freeMat(curJaco);
	freeMat(invOfJaco);
	freeMat(tsb);
	
	
	freeMat(workingMat);

}

int main(char* args) {
	//for the ABB IRB120
	float testHomeTFdata[] = {
	0, 0, 1, 374 ,
	0, 1, 0, 0,
	-1, 0,0, 630,
	0, 0, 0, 1
	};
	t_matrix testHomeTF =
	{ testHomeTFdata,4,4 };

	float testDesiredTFdata[] = {
	0.0125600000000000,0.811060000000000,- 0.584830000000000,- 179.686170000000,
- 0.388980000000000,	0.542770000000000,	0.744380000000000,	400,
0.921160000000000,	0.218140000000000,	0.322310000000000,	375.343850000000,
0,	0,	0,	1
	};
	t_matrix testDesiredTF =
	{ testDesiredTFdata,4,4 };

	printMatrix(&testHomeTF, "home transform mat");
	float SlistData[6][6] = {
		{-1.570796326794897,0,0,0,0,290},
		{0,0,1.570796326794897,0,- 270,0},
		{1.570796326794897,0,0,-70,0,0},
		{-1.570796326794897,0,0,0,0,302},
		{1.570796326794897,0,0,0,0,0},
		{0,0,0,0,0,72},
	};

	t_matrix Slist[6] = {
		{SlistData[0],6,1},
		{SlistData[1],6,1},
		{SlistData[2],6,1},
		{SlistData[3],6,1},
		{SlistData[4],6,1},
		{SlistData[5],6,1}
	};
	t_matrix initAngles = createMatrix(6, 1);
	initAngles.pData[0] = 1.81267624967643;
	initAngles.pData[1] = 1.55812741662800;
	initAngles.pData[2] = -1.52091948780048;
	initAngles.pData[3] = 2.30872895908779;
	initAngles.pData[4] = 0.554425086996533;
	initAngles.pData[5] = -5.61939843953905;

	t_matrix outputJointAngles = createMatrix(6, 1);
	IKinSpace(outputJointAngles.pData, &testDesiredTF, Slist, &testHomeTF, &initAngles, 0.01, 0.01);
		

/*
* expecting something like
	1.81267624967643
1.55812741662800
-1.52091948780048
2.30872895908779
0.554425086996533
-5.61939843953905
*/
	//printMatrix(&testMatB, "B");
	printMatrix(&outputJointAngles, "Result");
	freeMat(outputJointAngles);
	return(0);
}