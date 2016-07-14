var pages =  [];

var limit = 5;

var posts = $("#posts").children();
var c = 0;
var s = "";
posts.each(function(){
  s+= $(this).html();
  c++;
  if(c==limit){

    pages[pages.length] = s;
    s = "";
    c=0;

  }
});
if(c>0)pages[pages.length] = s;
var current_page = 0;
function setPage(p){


  $("#posts").fadeOut("slow",
  function(){
    $("#posts").html(pages[p]);
    current_page = p;
    $( "#next_page" ).prop( "disabled", !(current_page< pages.length-1) );
    $( "#prev_page" ).prop( "disabled", !(current_page>0) );
    $("#posts").fadeIn("slow",
      function(){
      }
    );
  });
}
setPage(0);
$("#next_page").click(function(){setPage(current_page+1);})
$("#prev_page").click(function(){setPage(current_page-1);})
