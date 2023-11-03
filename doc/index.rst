Biolearn
========

Biolearn enables easy and versatile analyses of biomarkers of aging data. It provides tools to easy load data from publicly available sources like the 
`Gene Expression Omnibus <https://www.ncbi.nlm.nih.gov/geo/>`_, `National Health and Nutrition Examimation Survey <https://www.cdc.gov/nchs/nhanes/index.htm>`_,
and the `Framingham Heart Study <https://www.framinghamheartstudy.org/>`_. Biolearns also contains reference implemenations for common aging clock such at the 
Horvath clock, DunedinPACE and many others that can easily be run in only a few lines of code.

Biolearn is scientifically and financially supported by `Biomarkers of Aging Consortium <https://www.agingconsortium.org/>`_ 
and `Methuselah Foundation <https://www.mfoundation.org/>`_ respectively 

.. grid::

    .. grid-item-card:: :fas:`rocket` Quickstart
        :link: quickstart
        :link-type: ref
        :columns: 12 12 4 4
        :class-card: sd-shadow-md
        :class-title: sd-text-primary
        :margin: 2 2 0 0

        Get started with Biolearn

    .. grid-item-card:: :fas:`th` Examples
        :link: auto_examples/index.html
        :link-type: url
        :columns: 12 12 4 4
        :class-card: sd-shadow-md
        :class-title: sd-text-primary
        :margin: 2 2 0 0

        Discover functionalities by reading examples

Featured examples
-----------------

.. grid::

  .. grid-item-card::
    :link: auto_examples/00_epigenetic_biomarkers/plot_epigenetic_clocks_on_geo.html
    :link-type: url
    :columns: 12 12 12 12
    :class-card: sd-shadow-sm
    :margin: 2 2 auto auto

    .. grid::
      :gutter: 3
      :margin: 0
      :padding: 0

      .. grid-item::
        :columns: 12 4 4 4

        .. image:: auto_examples/00_epigenetic_biomarkers/images/sphx_glr_plot_epigenetic_clocks_on_geo_001.png
          

      .. grid-item::
        :columns: 12 8 8 8

        .. div:: sd-font-weight-bold

          Demonstrate computation of several epigenetic clocks

        Show how the clocks compare with chronological age

  .. grid-item-card::
    :link: auto_examples/01_composite_biomarkers/plot_nhanes.html
    :link-type: url
    :columns: 12 12 12 12
    :class-card: sd-shadow-sm
    :margin: 2 2 auto auto

    .. grid::
      :gutter: 3
      :margin: 0
      :padding: 0

      .. grid-item::
        :columns: 12 4 4 4

        .. image:: auto_examples/01_composite_biomarkers/images/sphx_glr_plot_nhanes_002.png
          

      .. grid-item::
        :columns: 12 8 8 8

        .. div:: sd-font-weight-bold

          NHANES survival plotting with biological age

        Using NHANES data shows relationship between various data points and survival


.. toctree::
   :hidden:
   :includehidden:
   :titlesonly:

   quickstart.md
   clocks.rst
   data.rst
   auto_examples/index.rst
   modules/index.rst

.. toctree::
   :hidden:
   :caption: Development

   authors.rst
   GitHub Repository <https://github.com/bio-learn/biolearn>
