c = get_config()

c.NbConvertApp.notebooks = ['NRmin-report.ipynb']
c.NbConvertApp.export_format = 'latex'
c.NbConvertApp.postprocessor_class = 'PDF'

c.Exporter.template_file = 'NRmin-report'
