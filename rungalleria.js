Galleria.loadTheme('/galleria/themes/classic/galleria.classic.min.js');
Galleria.configure({
    showInfo: true,
    _toggleInfo: false
});
Galleria.ready(function() {
  var gallery = this;
  var fscr = this.$('images')
        .click(function() {
            gallery.toggleFullscreen();
        });
});
Galleria.run('.galleria');