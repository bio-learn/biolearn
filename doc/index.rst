Biolearn
========

Biolearn enables easy and versatile analyses of biomarkers of aging data. It provides tools to easily load data from publicly available sources like the 
`Gene Expression Omnibus <https://www.ncbi.nlm.nih.gov/geo/>`_, `National Health and Nutrition Examination Survey <https://www.cdc.gov/nchs/nhanes/index.htm>`_,
and the `Framingham Heart Study <https://www.framinghamheartstudy.org/>`_. Biolearn also contains reference implementations for common aging clocks such as the 
Horvath clock, DunedinPACE, and many others that can easily be run in only a few lines of code. You can read more about it in our `paper <https://www.biorxiv.org/content/10.1101/2023.12.02.569722v2>`_.

Biolearn is developed and supported by several organizations and individuals, especially `Biomarkers of Aging Consortium <https://www.agingconsortium.org/>`_,
`Methuselah Foundation <https://www.mfoundation.org/>`_, and `VOLO Foundation <https://www.VOLOfoundation.org/>`_.

We are hosting a 2024-2025 Challenge series on the Synapse platform, where participants will be asked to predict chronological age, 
mortality, and multi-morbidity, with total awards of $200k+. `Learn more at Synapse! <https://www.synapse.org/#!Synapse:syn52966292/wiki/624696/>`_.

If you use Biolearn in your research, please cite our `preprint <https://doi.org/10.1101/2023.12.02.569722>`_:

Ying, K., Paulson, S., Perez-Guevara, M., Emamifar, M., Martí nez, M. C., Kwon, D., Poganik, J. R., Moqri, M., & Gladyshev, V. N. (2023). Biolearn, an open-source library for biomarkers of aging. bioRxiv. https://doi.org/10.1101/2023.12.02.569722

.. grid::

    .. grid-item-card:: :fas:`rocket` Quickstart
        :link: quickstart
        :link-type: ref
        :columns: 12 12 4 4
        :class-card: sd-shadow-md
        :class-title: sd-text-primary
        :margin: 2 2 0 0

        Get started with Biolearn

    .. grid-item-card:: :fas:`road` Roadmap
        :link: roadmap.html
        :link-type: url
        :columns: 12 12 4 4
        :class-card: sd-shadow-md
        :class-title: sd-text-primary
        :margin: 2 2 0 0

        See where Biolearn is headed and how to contribute

    .. grid-item-card:: :fas:`trophy` Challenge
        :link: https://www.synapse.org/#!Synapse:syn52966292/
        :link-type: url
        :columns: 12 12 4 4
        :class-card: sd-shadow-md
        :class-title: sd-text-primary
        :margin: 2 2 0 0

        Join the Biomarkers of Aging Challenge!

Featured examples
-----------------

.. grid::

  .. grid-item-card::
    :link: auto_examples/00_omics_biomarkers/plot_epigenetic_clocks_on_geo.html
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

        .. image:: auto_examples/00_omics_biomarkers/images/sphx_glr_plot_epigenetic_clocks_on_geo_001.png
          

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
