<h1>imageTools</h1>

<h1><img src="logo_1.png" alt="imageTools Logo" style="width: 100px;"/> </h1>

<p>The <strong>imageTools</strong> package provides a suite of tools designed for the analysis and manipulation of raster images. This package includes functionalities for image classification, confusion matrix evaluation, data frame to image conversion, and more. Whether you're working with satellite imagery <strong>imageTools</strong> offers some utilities you need to streamline your workflow.</p>

<h2>Features</h2>
<ul>
    <li><strong>Land Cover Classification:</strong> Tools for LULC image classification using some statistical models (Maximum Likelihood Classifier).</li>
    <li><strong>Trajectory Classification:</strong> Tools for image trajectory classification using CMAP (Reis et al., 2020) approach (Compound Maximum a Posteriori).</li>
    <li><strong>Data Analysis:</strong> Functions created for data analysis such as boxplot creation, basic statistical report, class separability analysis by stochastic distance, and dendrogram creation, among others.</li>
    <li><strong>Confusion Matrix Evaluation:</strong> Functions to evaluate classification performance through confusion matrices.</li>
    <li><strong>Image Manipulation:</strong> Functions to convert images to data frames and vice versa.</li>
    <li><strong>Trejactory Analysis:</strong> Tools to calculate and visualize impossible transitions in time series of raster classifications.</li>
</ul>

<h2>Installation</h2>
<p>To install the latest version of the <strong>imageTools</strong> package from GitHub, use the following commands in R:</p>

<pre><code>install.packages("devtools")
devtools::install_github("yourusername/imageTools")</code></pre>

<h2>Usage</h2>

<h3>Data Frame to Image Conversion</h3>
<p>Convert a data frame to an image array, which is useful for manipulating raster data as data frames and vice versa:</p>

<pre><code># Load the necessary libraries
library(raster)
# Create a data frame with random values
df_img <- data.frame(Band1 = runif(100), Band2 = runif(100), Band3 = runif(100))

# Create a RasterStack with random values
template_raster <- stack(replicate(3, raster(matrix(runif(100), 10, 10))))

# Define band names
band_names <- c("Band1", "Band2", "Band3")

# Convert image to data frame
df_img <- img2df(template_raster, band_names)

# Print the resulting data frame
print(df_img)</code></pre>

<h2>Documentation</h2>
<p>Detailed documentation for each function is available within the package. You can access the documentation using the <code>help</code> function in R. For example:</p>

<pre><code>help(CMAP_classifier)
help(evaluate_classification)</code></pre>

<h2>Contributing</h2>
<p>Contributions to the <strong>imageTools</strong> package are welcome. If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request on GitHub. You can also contact me by email (vinicius.dlucas@gmail.com) </p>

<h2>License</h2>
<p>This package is licensed under the MIT License. See the <code>LICENSE</code> file for more details.</p>

</body>
</html>
