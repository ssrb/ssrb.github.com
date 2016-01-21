---
layout: default
title:
tagline:
---

<div class="row-fluid" id="page-index">
	<div class="span8">
		<div class="posts">
		{% for post in site.posts %}
 		{% if post.staging != true%}
		<article class="post">
		    <header>
		      <h2><a href="{{ BASE_PATH }}{{ post.url }}">{{ post.title }}</a></h2>
		      <h3>{{ post.date | date_to_string }}</h3>
		    </header>
		    {{ post.content | split: '<!-- more -->' | first }}
		    <a href="{{ BASE_PATH }}{{ post.url }}">more ...</a>
		</article>
		<hr/>
 		{% endif %}
		{% endfor %}
		</div>
	</div>

	<div class="span4 sidebar">
	  <section class="about">
	  <!--img class="icon" src="http://en.gravatar.com/userimage/18679074/9775b23e0c4499773692e8b2f8060b93.png?size=100"-->
	  <h2>
	    cat ~/.plan
	  </h2>
	  <p>
	  	I'm a middleware developer based in Christchurch, New Zealand.
	  </p>
	  <p>
	  	I worked for different industries and several research institutes. 
	  	If you're interested in the details, feel free to take a look at my <a href="/cv_sbigot.pdf">CV</a>.
	  </p>
	  <p>
	  	These days, I'm enjoying my time with computer engineering and electronics, computer science and mathematics, hacking and making all sort of things.
	  </p>
	  </section>
	</div>
</div>

