package edu.tum.cs.influence.meso;

import java.util.Arrays;
import java.util.List;
import java.util.Queue;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.OptimizationRequest;
import com.joptimizer.optimizers.OptimizationResponse;

import edu.tum.cs.math.dist.DiscreteDistribution;

public class PredictiveModel {

	private static class SumOfJensenShannonDivergencesSimple implements ConvexMultivariateRealFunction {
		private final double[][][] predictors;	// instance x dimension x variable
		private final double[][] response;	// instance x dimension

		public SumOfJensenShannonDivergencesSimple(double[][][] predictors, double[][] response) {
			this.predictors = predictors;
			this.response = response;
		}

		@Override
		public int getDim() {
			return predictors[0][0].length;
		}

		@Override
		public double value(double[] c) {
			double sumJsd = 0.0;
			for (int i = 0; i < predictors.length; i++) {
				for (int j = 0; j < predictors[0].length; j++) {
					double sumComp = 0.0;
					for (int k = 0; k < predictors[0][0].length; k++)
						sumComp += c[k] * predictors[i][j][k];
					sumJsd += ((DiscreteDistribution.LOG2 - Math.log(sumComp + response[i][j])) *
							(sumComp + response[i][j])) + (Math.log(sumComp) * sumComp);
				}
			}
			return sumJsd;
		}

		@Override
		public double[] gradient(double[] c) {
			double[] grad = new double[c.length];
			for (int i = 0; i < predictors.length; i++) {
				for (int j = 0; j < predictors[0].length; j++) {
					double sumComp = 0.0;
					for (int k = 0; k < predictors[0][0].length; k++)
						sumComp += c[k] * predictors[i][j][k];
					double v = DiscreteDistribution.LOG2 - Math.log(sumComp + response[i][j]) + Math.log(sumComp);
					for (int k = 0; k < predictors[0][0].length; k++)
						grad[k] += v * predictors[i][j][k];
				}
			}
			return grad;
		}

		@Override
		public double[][] hessian(double[] c) {
			double[][] h = new double[c.length][c.length];
			for (int i = 0; i < predictors.length; i++) {
				for (int j = 0; j < predictors[0].length; j++) {
					double sumComp = 0.0;
					for (int k = 0; k < predictors[0][0].length; k++)
						sumComp += c[k] * predictors[i][j][k];
					double v = response[i][j] / ((sumComp + response[i][j]) * sumComp);
					for (int row = 0; row < predictors[0][0].length; row++) {
						for (int col = 0; col <= row; col++) {
							double w = v * predictors[i][j][row] * predictors[i][j][col];
							h[row][col] += w;
							if (col != row)
								h[col][row] += w;
						}
					}
				}
			}
			return h;
		}
	}

	private static class SumOfJensenShannonDivergences implements ConvexMultivariateRealFunction {
		private final double[][][] predictors;	// instance x dimension x variable
		private final double[][] response;	// instance x dimension
		private final double alphaL1, alphaL2;
		private final int numCoeffParam;

		public SumOfJensenShannonDivergences(double[][][] predictors, double[][] response, double alphaL1,
				double alphaL2) {
			this.predictors = predictors;
			this.response = response;
			this.alphaL1 = alphaL1;
			this.alphaL2 = alphaL2;
			numCoeffParam = predictors[0][0].length;
		}

		@Override
		public int getDim() {
			return numCoeffParam + predictors[0].length;
		}

		@Override
		public double value(double[] c) {
			double sumJsd = 0.0;
			for (int i = 0; i < predictors.length; i++) {
				double norm1 = 0.0, norm2 = 0.0;
				for (int j = 0; j < predictors[0].length; j++) {
					double sumComp = c[numCoeffParam + j];
					for (int k = 0; k < predictors[0][0].length; k++)
						sumComp += c[k] * predictors[i][j][k];
					norm1 += sumComp;
					norm2 += sumComp * sumComp;
					sumJsd += ((DiscreteDistribution.LOG2 - Math.log(sumComp + response[i][j])) *
							(sumComp + response[i][j])) + (Math.log(sumComp) * sumComp);
				}

				// L1 regularization
				// This simplification relies on the IP solver to never evaluate the function or its derivatives for
				// values that violate the constraints, specifically the non-negativity.
				sumJsd += alphaL1 * norm1;

				// L2 regularization
				sumJsd += alphaL2 * Math.sqrt(norm2);
			}
			return sumJsd;
		}

		@Override
		public double[] gradient(double[] c) {
			double[] grad = new double[c.length];
			double[] sumComp = new double[predictors[0].length];
			double[] norm1 = new double[predictors[0][0].length];
			double[] norm2Num = new double[predictors[0][0].length];
			for (int i = 0; i < predictors.length; i++) {
				Arrays.fill(sumComp, 0.0);
				Arrays.fill(norm1, 0.0);
				Arrays.fill(norm2Num, 0.0);
				double norm2Denom = 0.0;

				for (int j = 0; j < predictors[0].length; j++) {
					sumComp[j] = c[numCoeffParam + j];
					for (int k = 0; k < predictors[0][0].length; k++)
						sumComp[j] += c[k] * predictors[i][j][k];
					double v = DiscreteDistribution.LOG2 - Math.log(sumComp[j] + response[i][j]) +
							Math.log(sumComp[j]);

					// coefficient parameters
					for (int k = 0; k < predictors[0][0].length; k++) {
						norm1[k] += predictors[i][j][k];
						norm2Num[k] += sumComp[j] * predictors[i][j][k];
						grad[k] += v * predictors[i][j][k];
					}
					norm2Denom += sumComp[j] * sumComp[j];

					// non-coefficient parameters
					grad[numCoeffParam + j] += v;
				}

				// L1/L2 regularization
				norm2Denom = Math.sqrt(norm2Denom);
				int idx = 0;
				for (int k = 0; k < predictors[0][0].length; k++) {
					grad[idx] += alphaL1 * norm1[k];
					grad[idx] += alphaL2 * (norm2Num[k] / norm2Denom);
					idx++;
				}
				for (int k = 0; k < predictors[0].length; k++) {
					grad[idx] += alphaL1;
					grad[idx] += alphaL2 * (sumComp[k] / norm2Denom);
					idx++;
				}
			}
			return grad;
		}

		@Override
		public double[][] hessian(double[] c) {
			double[][] h = new double[c.length][c.length];
			double[] sumComp = new double[predictors[0].length];
			double[] normSum1 = new double[c.length];
			double[][] normSum12 = new double[c.length][c.length];
			for (int i = 0; i < predictors.length; i++) {
				double normSum = 0.0;
				Arrays.fill(normSum1, 0.0);
				for (int j = 0; j < c.length; j++)
					Arrays.fill(normSum12[j], 0.0);

				for (int j = 0; j < predictors[0].length; j++) {
					sumComp[j] = c[numCoeffParam + j];
					for (int k = 0; k < predictors[0][0].length; k++)
						sumComp[j] += c[k] * predictors[i][j][k];
					double v = response[i][j] / ((sumComp[j] + response[i][j]) * sumComp[j]);

					// top left (dense, symmetric): coefficient x coefficient
					for (int row = 0; row < predictors[0][0].length; row++) {
						normSum1[row] += sumComp[j] * predictors[i][j][row];
						for (int col = 0; col <= row; col++) {
							double w = predictors[i][j][row] * predictors[i][j][col];
							normSum12[row][col] += w;
							w *= v;
							h[row][col] += w;
							if (col != row)
								h[col][row] += w;
						}
					}
					normSum += sumComp[j] * sumComp[j];

					// bottom left (dense): non-coefficient x coefficient
					for (int col = 0; col < predictors[0][0].length; col++) {
						h[numCoeffParam + j][col] += v * predictors[i][j][col];
						h[col][numCoeffParam + j] += v * predictors[i][j][col];
					}

					// bottom right (diagonal): non-coefficient x non-coefficient
					h[numCoeffParam + j][numCoeffParam + j] += v;
				}

				// L1 regularization term is 0
				// L2 regularization top left
				double normSumSqrt = Math.sqrt(normSum);
				for (int row = 0; row < predictors[0][0].length; row++) {
					for (int col = 0; col <= row; col++) {
						double w = alphaL2 * ((normSum12[row][col] * (normSumSqrt / normSum)) -
								((normSum1[row] * normSum1[col]) / (normSumSqrt * normSum)));
						h[row][col] += w;
						if (col != row)
							h[col][row] += w;
					}
				}

				// L2 regularization bottom
				for (int row = predictors[0][0].length; row < h.length; row++) {
					// bottom left
					for (int col = 0; col < predictors[0][0].length; col++) {
						double w = alphaL2 * ((predictors[i][row - predictors[0][0].length][col] *
									(normSumSqrt / normSum)) -
								((sumComp[row - predictors[0][0].length] * normSum1[col]) / (normSumSqrt * normSum)));
						h[row][col] += w;
						h[col][row] += w;
					}

					// bottom right (excluding the diagonal)
					for (int col = predictors[0][0].length; col < row; col++) {
						double w = alphaL2 * (-(sumComp[row - predictors[0][0].length] *
									sumComp[col - predictors[0][0].length]) / (normSumSqrt * normSum));
						h[row][col] += w;
						h[col][row] += w;
					}

					// bottom diagonal
					h[row][row] += alphaL2 * ((normSum - (sumComp[row - predictors[0][0].length] *
							sumComp[row - predictors[0][0].length])) / (normSumSqrt * normSum));
				}
			}
			return h;
		}
	}

	private final int numTopics;
	private final int numCoefficients;
	private final double alphaL1, alphaL2;
	private double[] coefficients;

	private PredictiveModel(int numTopics, int numCoefficients, double alphaL1, double alphaL2, double[] coefficients) {
		if (coefficients.length != numCoefficients)
			throw new IllegalArgumentException("invalid number of coefficients");

		this.numTopics = numTopics;
		this.numCoefficients = numCoefficients;
		this.alphaL1 = alphaL1;
		this.alphaL2 = alphaL2;
		this.coefficients = coefficients;
	}

	public PredictiveModel(int numTopics, int numCoefficients, double alphaL1, double alphaL2) {
		this(numTopics, numCoefficients + numTopics, alphaL1, alphaL2, new double[numCoefficients + numTopics]);
		int numTrueCoeff = coefficients.length - numTopics;
		for (int i = 0; i < numTrueCoeff; i++)
			coefficients[i] = 1.0 / numTrueCoeff;
	}

	public double[] getCoefficients() {
		return coefficients;
	}

	private double[] minimizeFunction(ConvexMultivariateRealFunction f0) {
		int numParams = f0.getDim();

		// equality constraints
		double[][] a = new double[1][numParams];
		Arrays.fill(a[0], 1.0);
		double[] b = new double[] { 1.0 };

		// inequality constraints
		ConvexMultivariateRealFunction[] fi = new ConvexMultivariateRealFunction[numParams];
		for (int i = 0; i < numParams; i++) {
			double[] q = new double[numParams];
			q[i] = -1.0;
			fi[i] = new LinearMultivariateRealFunction(q, 0);
		}

		// perform optimization
		OptimizationRequest req = new OptimizationRequest();
		req.setMaxIteration(10 * req.getMaxIteration());
		req.setF0(f0);
		req.setFi(fi);
		req.setA(a);
		req.setB(b);
		double[] initialPoint = new double[numParams];
		Arrays.fill(initialPoint, 1.0 / numParams);
		req.setInitialPoint(initialPoint);
		JOptimizer opt = new JOptimizer();
		opt.setOptimizationRequest(req);
		try {
			opt.optimize();
		} catch (Exception ex) {
			throw new RuntimeException("optimizer threw exception", ex);
		}
		OptimizationResponse res = opt.getOptimizationResponse();
		if (res.getReturnCode() == OptimizationResponse.FAILED)
			throw new RuntimeException("optimization failed");
		return res.getSolution();
	}

	public SummaryStatistics fitParameters(List<InfluencedUser> trainUsers) {
		if (trainUsers.isEmpty())
			throw new RuntimeException("empty training set");

		double[][][] x = new double[trainUsers.size()][][];
		double[][] r = new double[trainUsers.size()][];
		int idx = 0;
		for (InfluencedUser user : trainUsers) {
			x[idx] = user.getThetaMatrix();
			r[idx] = user.getThetaToPredict();
			idx++;
		}

		// minimize sum of JSDs
		SumOfJensenShannonDivergences f0 = new SumOfJensenShannonDivergences(x, r, alphaL1, alphaL2);
		coefficients = minimizeFunction(f0);

		// compute training error
		SummaryStatistics errorStats = new SummaryStatistics();
		for (InfluencedUser user : trainUsers)
			errorStats.addValue(DiscreteDistribution.distJS2(predictTheta(user), user.getThetaToPredict()));
		return errorStats;
	}

	public SummaryStatistics computePerUserError(List<InfluencedUser> trainUsers) {
		if (trainUsers.isEmpty())
			throw new RuntimeException("empty training set");

		SummaryStatistics errorStats = new SummaryStatistics();
		double[][][] x = new double[1][][];
		double[][] r = new double[1][];
		for (InfluencedUser user : trainUsers) {
			x[0] = user.getThetaMatrix();
			r[0] = user.getThetaToPredict();
			// using optimization without additional component and without regularization
			SumOfJensenShannonDivergencesSimple f0 = new SumOfJensenShannonDivergencesSimple(x, r);
			double[] coefficients = minimizeFunction(f0);
			errorStats.addValue(DiscreteDistribution.distJS2(predictTheta(user, coefficients, true),
					user.getThetaToPredict()));
		}
		return errorStats;
	}

	private double[] predictTheta(InfluencedUser user, double[] coefficients, boolean isSimple) {
		double[] predictedTheta = new double[numTopics];
		double sum = 0.0;

		double[][] currentThetas = user.getThetaMatrix();
		for (int i = 0; i < currentThetas.length; i++) {
			for (int j = 0; j < currentThetas[i].length; j++)
				predictedTheta[i] += coefficients[j] * currentThetas[i][j];
			if (!isSimple)
				predictedTheta[i] += coefficients[currentThetas[i].length + i];
			sum += predictedTheta[i];
		}

		if (Math.abs(1.0 - sum) > InfluencedUser.EPSILON)
			throw new RuntimeException("sum of predicted edge theta != 1: " + sum);
		return predictedTheta;
	}

	public double[] predictTheta(InfluencedUser user) {
		return predictTheta(user, coefficients, false);
	}

	private static final int maxCoefficients = 14 + 150;

	public static String getCsvHeader() {
		StringBuilder sb = new StringBuilder();
		sb.append("numTopics;numCoefficients;alphaL1;alphaL2");
		for (int i = 0; i < maxCoefficients; i++)
			sb.append(';').append('c').append(i + 1);
		return sb.toString();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(numTopics).append(';').append(numCoefficients).append(';');
		sb.append(alphaL1).append(';').append(alphaL2);
		for (int i = 0; i < maxCoefficients; i++) {
			sb.append(';');
			if (i < numCoefficients)
				sb.append(coefficients[i]);
			else
				sb.append('0');
		}
		return sb.toString();
	}

	public static PredictiveModel readCsv(Queue<String> parts) {
		int numTopics = Integer.parseInt(parts.poll());
		int numCoefficients = Integer.parseInt(parts.poll());
		double alphaL1 = Double.parseDouble(parts.poll());
		double alphaL2 = Double.parseDouble(parts.poll());
		double[] coefficients = new double[numCoefficients];
		for (int i = 0; i < maxCoefficients; i++) {
			String coeffStr = parts.poll();
			if (i < numCoefficients)
				coefficients[i] = Double.parseDouble(coeffStr);
		}
		return new PredictiveModel(numTopics, numCoefficients, alphaL1, alphaL2, coefficients);
	}

}
