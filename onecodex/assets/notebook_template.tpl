{%- extends 'full.tpl' -%}

{%- block html_head -%}
<meta charset="utf-8" />
<div id="reportdate">{{ resources['metadata']['date'] }}</div>
{% set nb_title = nb.metadata.get('title', '') or resources['metadata']['name'] %}
<title>{{nb_title}}</title>
<link rel="stylesheet" href="onecodex/assets/notebook_template.css">
{%- endblock html_head -%}

{% block output_area_prompt %}
{% endblock output_area_prompt %}

{% block input_group %}
<div class="input_hidden">
  {% if cell.metadata.show_input == True %}
    {{ super() }}
  {% endif %}
</div>
{% endblock input_group %}

{% block output_group %}
<div class="{{cell.metadata.classes}}">
  {{ super() }}
</div>
{% endblock %}
