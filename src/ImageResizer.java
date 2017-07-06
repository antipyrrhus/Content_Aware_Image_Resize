/** Class: ImageResizer.java
 *  @author Yury Park
 *
 *  Seam-carving is a content-aware image resizing technique where the image is reduced in size by one pixel of height (or width)
 *  at a time. A vertical seam in an image is a path of pixels connected from the top to the bottom with one pixel in each row;
 *  a horizontal seam is a path of pixels connected from the left to the right with one pixel in each column.
 *
 *  1. Energy calculation. The first step is to calculate the energy of a pixel, which is a measure of its importance --
 *  the higher the energy, the less likely that the pixel will be included as part of a seam (as we'll see in the next step).
 *  The energy is high (white) for pixels in the image where there is a rapid color gradient. The seam-carving technique
 *  avoids removing such high-energy pixels.
 *
 *  2. Seam identification. The next step is to find a vertical seam of minimum total energy. This is similar to the
 *  classic shortest path problem in an edge-weighted digraph, but there are three important differences:
 *  - The weights are on the vertices instead of the edges.
 *  - The goal is to find the shortest path from any of the W pixels in the top row to any of the W pixels in the bottom row.
 *  - The digraph is acyclic, where there is a downward edge from pixel (x, y) to pixels (x - 1, y + 1), (x, y + 1), and
 *    (x + 1, y + 1). assuming that the coordinates are in the prescribed ranges. Seams cannot wrap around the image
 *    (e.g., a vertical seam cannot cross from the leftmost column of the image to the rightmost column).
 *  Finding a horizontal seam is analogous.
 *
 *  3. Seam removal. The final step is remove from the image all of the pixels along the vertical or horizontal seam.
 *
 */
import java.util.*;
import java.awt.Color;

public class ImageResizer {
    private Picture pic;
    private final int BORDER_PIXEL_ENERGY = 1000;   //all pixels at the outer border are set to 1000 by default
    private int[][] picColors;      //color of each pixel at (x, y). This is effectively a backup in case this.pic is edited.
    private boolean isTransposed;   //keeps track of whether the picture is transposed or not.

    /**
     * 1-arg constructor.
     * @param picture
     */
    public ImageResizer(Picture picture)       // create an image resizer object based on the given picture
    {
        checkIfNull(picture);   //throw exception if the parameter is null
        this.pic = picture;
        this.isTransposed = false;   //picture is originally not transposed

        /* Save each pixel color in 2D array. */
        this.picColors = new int[width()][height()];
        for (int x = 0; x < picture.width(); x++) {
            for (int y = 0; y < picture.height(); y++) {
                this.picColors[x][y] = pic.get(x, y).getRGB();
            }
        }
    }

    /**
     * Method: picture
     * @return the picture
     */
    public Picture picture()                          // current picture
    {
        //If picture is transposed, transpose it back to the original orientation
    	//This also resets the picColors array, which will be relevant below.
        if (this.isTransposed) {
            this.transposePic();
            this.isTransposed = false;
        }

        /* Note: What if the end user makes edits to the picture thru some 3rd-party app
         * after having constructed this class? We want to make sure that any outside
         * edits like that shouldn't affect the picture that we return.
         * So, as a safeguard, we will construct a new picture from
         * the array of colors that we have been keeping as a backup separately, and return that. */
        this.pic = new Picture(width(), height());
        for (int x = 0; x < width(); x++) {
            for (int y = 0; y < height(); y++) {
                this.pic.set(x, y, new Color(this.picColors[x][y]));
            }
        }
        return pic;
    }

    /**
     * Method: width
     * @return picture's width
     */
    public int width()                            // width of current picture
    {
        return pic.width();
    }

    /**
     * Method: height
     * @return picture's height
     */
    public int height()                           // height of current picture
    {
        return pic.height();
    }

    /**
     * Method: energy
     * @param x column index
     * @param y row index
     * @return the energy of pixel at the given column and row.
     */
    public double energy(int x, int y)               // energy of pixel at column x and row y
    {
    	//Check for illegal parameters
        if (x < 0 || y < 0 || x > width() - 1 || y > height() - 1) {
            throw new java.lang.IndexOutOfBoundsException();
        }

        /* Base case:
         * Energy of pixel at the border is set to 1000 by default so that it is strictly
         * larger than any of the inner pixels' energy.
         */
        if (x == 0 || y == 0 || x == width() - 1 || y == height() - 1) {
            return BORDER_PIXEL_ENERGY;
        }

        /* Custom method calls to compute energy at pixels that are inside the border */
        return Math.sqrt(deltaXSquared(x,y) + deltaYSquared(x,y));
    }

    /**
     * Method: deltaXSquared
     * @param x column index
     * @param y row index
     * @return double value after performing the calculations for pixel's energy at the given (x,y) location
     *         in terms of the gradient difference of its x+1 and x-1 neighboring pixels.
     */
    private double deltaXSquared(int x, int y) {
        int rx, gx, bx;
        Color pixel_at_x_plus1_comma_y  = pic.get(x+1, y); //java.awt.Color
        Color pixel_at_x_minus1_comma_y = pic.get(x-1, y);

        rx = pixel_at_x_plus1_comma_y.getRed()   - pixel_at_x_minus1_comma_y.getRed();
        gx = pixel_at_x_plus1_comma_y.getGreen() - pixel_at_x_minus1_comma_y.getGreen();
        bx = pixel_at_x_plus1_comma_y.getBlue()  - pixel_at_x_minus1_comma_y.getBlue();
        return (rx * rx) + (gx * gx) + (bx * bx);
    }

    /**
     * Method: deltaYSquared
     * @param x column index
     * @param y row index
     * @return double value. Calculation is similar to that in deltaXSquared except this time it's in terms of
     *         the y+1 and y-1 neighboring pixels.
     */
    private double deltaYSquared(int x, int y) {
        int ry, gy, by;
        java.awt.Color pixel_at_x_comma_y_plus1  = pic.get(x, y+1);
        java.awt.Color pixel_at_x_comma_y_minus1 = pic.get(x, y-1);

        ry = pixel_at_x_comma_y_plus1.getRed()   - pixel_at_x_comma_y_minus1.getRed();
        gy = pixel_at_x_comma_y_plus1.getGreen() - pixel_at_x_comma_y_minus1.getGreen();
        by = pixel_at_x_comma_y_plus1.getBlue()  - pixel_at_x_comma_y_minus1.getBlue();
        return (ry * ry) + (gy * gy) + (by * by);
    }

    /**
     * Method: findHorizontalSeam
     * @return a minimum-energy horizontal seam to be removed
     */
    public int[] findHorizontalSeam()   // sequence of indices for horizontal seam
    {
        return findSeam(true);  //custom method call
    }

    /**
     * Method: findVerticalSeam
     * @return a minimum-energy vertical seam to be removed
     */
    public int[] findVerticalSeam()     // sequence of indices for vertical seam
    {
        return findSeam(false); //custom method call
    }

    /**
     * Method: transposePic
     * Transpose the picture.
     */
    private void transposePic() {
        Picture transposedPic = new Picture(pic.height(), pic.width()); //initialize transposed version of the picture
        int[][] transposedPicColors = new int[height()][width()];       //initialize the color array of the transposed pic as a backup
        for (int x = 0; x < this.pic.width(); x++) {
            for (int y = 0; y < this.pic.height(); y++) {
                transposedPic.set(y, x, this.pic.get(x, y));
                transposedPicColors[y][x] = this.picColors[x][y];
            }
        }

        this.pic = transposedPic;
        this.picColors = transposedPicColors;
    }

    /**
     * Method: findSeam. Finds either a horizontal or vertical seam to be removed.
     * @param mustTransposePic whether the picture should be transposed from its original form.
     *        This effectively determines whether we're finding a vertical or horizontal seam,
     *        without having to change our method in any other way.
     * @return an int[] array consisting of the x- or y- locations of pixels consisting a seam.
     */
    private int[] findSeam(boolean mustTransposePic) {

        if (mustTransposePic != this.isTransposed) {
            transposePic();
            this.isTransposed = !this.isTransposed;
        }

        //Implement a shortest-path (minimum-energy) algorithm similar to an acyclic topological algorithm
        /* We'll initialize a double[] array that keeps track of the distance to a pixel (aka "vertex")
         * We will represent the vertex location (x,y) as a single integer instead since we're using 1D array here,
         * which is why we designate a space of width () * height() */
        double[] distTo = new double[width() * height()];
        Arrays.fill(distTo, Double.POSITIVE_INFINITY);  //Initialize all pixels as having infinite distance
        for (int i = 0; i < pic.width(); i++)
            distTo[i] = this.BORDER_PIXEL_ENERGY;   //Now RE-Initialize the first row of pixels as having distance 1000

        //Another array to keep track of parent vertices (for constructing optimal shortest path)
        int[] parentOf = new int[width() * height()];
        for (int i = 0; i < parentOf.length; i++)
            parentOf[i] = i;    //initialize every pixel vertex as its own parent.

        /* Since we stipulate that every pixel in a picture points downward in SW, S and SE directions
         * this picture, when thought of as a graph, is already topologically sorted!
         * So we just have to go thru every vertex from top row to bottom row - 1 (we don't need to
         * go thru the bottom row itself because the bottom row pixels are the destinations),
         * and "relax" each vertex, which means updating the min. total distance to each vertex
         * considering the distance travelled thus far plus the energy cost. */
        for (int y = 0; y < pic.height() - 1; y++) {
            for (int x = 0; x < pic.width(); x++) {
                relax(x, y, distTo, parentOf); //custom method
            }
        }

        /* Once all pixel vertices are "relaxed", every pixel will now have
         * a minimum distTo[] value. We want to find the minimum TOTAL distance (minimum energy)
         * path from the top row to the bottom row. So we'll start
         * at the last row, see which vertices in the last row have the smallest distTo[]
         * value, and then backtrack our way to the starting vertex at the top row. */
        int lastY = height() - 1;
        double minEnergy = Double.POSITIVE_INFINITY;    //We want to see which vertex has the min. energy.
        int minEnergyXY = -1;                           //The x,y location of the vertex sith the smallest distTo[]

        //Go thru every vertex in the last row and find the vertex with the min. energy.
        for (int x = 0; x < width(); x++) {
            //Represent the (x,y) location of the current vertex as a single integer
            //(as you would for a 1-D array representation of a 2-D array)
            int thisXY = lastY * width() + x;
            if (distTo[thisXY] < minEnergy) {
                minEnergy = distTo[thisXY];
                minEnergyXY = thisXY;
            }
        }

        //Now that we found the min. energy vertex in the last row,
        //we'll use the parentOf[] array to backtrack our way to the top row vertex,
        //thereby creating a path that represents a min-energy seam.
        Stack<Integer> st = new Stack<>();
        st.push(minEnergyXY % width()); //this % operation decodes the 1-D array representation of the (x,y) location and just computes the x-location of the vertex
        int currXY = minEnergyXY;

        //Build the min-energy seam path by backtracking.
        while (st.size() < height()) {
            st.push(parentOf[currXY] % width());
            currXY = parentOf[currXY];
        }

        //Now that we've built the min-energy path backwards, we want to pop it back out
        //so that the path will be in order.
        int index = 0;
        int[] ret = new int[height()];
        while (!st.isEmpty()) {
            ret[index++] = st.pop();
        }

        //Make sure the picture is reverted to its original dimension if necessary
        if (this.isTransposed) {
            transposePic();
            this.isTransposed = false;
        }
        return ret;
    }

    /**
     * Method: relax. Invoked by findSeam().
     *         Updates the min. total distance to the neighbors of the vertex at (x,y)
     *         considering the distance travelled thus far plus the energy cost
     * @param x
     * @param y
     * @param distTo
     * @param parentOf
     */
    private void relax(int x, int y, double[] distTo, int[] parentOf) {
        int thisXY = y * width() + x;   //1-D array (single integer) representation of (x, y)

        //For every pixel, we must check its SW, S and SE neighbor.
        int nextY = y + 1;
        int prevX = x - 1;
        int thisX = x;
        int nextX = x + 1;

        int nextYPrevX = nextY * width() + prevX; //Single-integer representation of (x-1, y+1)
        int nextYThisX = nextY * width() + thisX; //Single-integer representation of (x,   y+1)
        int nextYNextX = nextY * width() + nextX; //Single-integer representation of (x+1, y+1)

        //Relax the SW neighbor vertex
        if (prevX >= 0) {
            relax(thisXY, nextY, prevX, nextYPrevX, distTo, parentOf); //helper method
        }

        //Relax the S neighbor vertex
        relax(thisXY, nextY, thisX, nextYThisX, distTo, parentOf); //helper method

        //Relax the SE neighbor vertex
        if (nextX < width()) {
            relax(thisXY, nextY, nextX, nextYNextX, distTo, parentOf); //helper method
        }
    }

    /**
     * Method: relax. Helper method invoked by relax() method above.
     * @param thisXY
     * @param otherY
     * @param otherX
     * @param otherXY
     * @param distTo
     * @param parentOf
     */
    private void relax(int thisXY, int otherY, int otherX, int otherXY, double[] distTo, int[] parentOf) {
        //Update the distTo[] and parentOf[] value for the vertex located at otherXY if appropriate
        double dist = distTo[thisXY] + this.energy(otherX, otherY);
        if (dist < distTo[otherXY]) {
            distTo[otherXY] = dist;
            parentOf[otherXY] = thisXY;
        }
    }

    /* Throw a java.lang.IllegalArgumentException if removeVerticalSeam() or removeHorizontalSeam() is
     * called with an array of the wrong length or if the array is not a valid seam
     * (i.e., either an entry is outside its prescribed range or two adjacent entries differ by more than 1).
     *
     * Throw a java.lang.IllegalArgumentException if removeVerticalSeam() is called when the width of
     * the picture is less than or equal to 1 or if removeHorizontalSeam() is called when the height
     * of the picture is less than or equal to 1. */
    /**
     * Method: removeHorizontalSeam. Invoked by ResizeDemo class.
     * @param seam
     */
    public    void removeHorizontalSeam(int[] seam)   // remove horizontal seam from current picture
    {
        checkIfNull(seam);
        //If the picture has a height of just 1 pixel, it doesn't make sense to remove a horizontal seam.
        if (height() <= 1) throw new java.lang.IllegalArgumentException();
        checkIfValid(seam, true);

        //Construct a new picture that will be the result of having horizontal seam removed
        Picture newPic = new Picture(width(), height() - 1);
        picColors = new int[width()][height() - 1]; //We'll also update the picColors array[][]

        //We'll copy the contents of the original picture to the newPic, except we won't copy
        //the contents of the horizontal seam. So we'll keep a separate track of the indices
        //of the newPic (newX and newY)
        int newX = 0;
        int newY = 0;
        for (int x = 0; x < width(); x++) {
            for (int y = 0; y < height(); y++) {
                if (y == seam[x]) continue; //If (x,y) is the location of the removed seam, don't copy it to the newPic
                newPic.set(newX, newY, pic.get(x, y));
                this.picColors[newX][newY] = newPic.get(newX, newY).getRGB();
                newY++;
            }
            newY = 0;
            newX++;
        }
        pic = newPic;
    }

    /**
     * Method: removeVerticalSeam. Invoked by ResizeDemo class.
     * @param seam
     */
    public    void removeVerticalSeam(int[] seam)     // remove vertical seam from current picture
    {
        checkIfNull(seam);
        if (width() <= 1) throw new java.lang.IllegalArgumentException();
        checkIfValid(seam, false);

        //Same logic as the removeHorizontalSeam() method.
        Picture newPic = new Picture(width() - 1, height());
        picColors = new int[width() - 1][height()];
        int newX = 0;
        int newY = 0;
        for (int y = 0; y < height(); y++) {
            for (int x = 0; x < width(); x++) {
                if (x == seam[y]) continue;
                newPic.set(newX, newY, pic.get(x, y));
                this.picColors[newX][newY] = newPic.get(newX, newY).getRGB();
                newX++;
            }
            newX = 0;
            newY++;
        }
        pic = newPic;
    }

    /**
     * Method: checkIfValid. Checks if the seam has correct length, within the correct bounds,
     *         and neighboring pixels of the seam differ by at most 1 pixel.
     * @param seam
     * @param isHorizontal whether the seam is horizontal or vertical.
     */
    private void checkIfValid(int[] seam, boolean isHorizontal) {
        //Check whether seam has correct length
        if (seam.length != (isHorizontal ? width() : height()))
            throw new java.lang.IllegalArgumentException();

        for (int i = 0; i < seam.length; i++) {
            //Check if the contents of seam are within correct bounds.
            if (seam[i] < 0 || seam[i] > (isHorizontal? height() - 1 : width() - 1))
                throw new java.lang.IllegalArgumentException();

            //Check if the contents of seam differ by at most 1 pixel from one to the next.
            //If not, not a valid seam
            if (i == 0) continue;
            if (Math.abs(seam[i] - seam[i-1]) > 1)
                throw new java.lang.IllegalArgumentException();
        }
    }

    /**
     * Method: checkIfNull
     * If the parameter is null, throw exception.
     *
     * @param object any object
     */
    private void checkIfNull(Object object) {
        if (object == null) throw new java.lang.NullPointerException();
    }

    /* NO MAIN CLASS. RUN ResizeDemo.java. */
 }