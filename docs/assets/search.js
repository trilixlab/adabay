/* Tiny client-side search for the adabay docs site.
 * Looks up the index JSON relative to the data-base attribute on the
 * search input, then filters by substring match against title + section + kw.
 */
(function () {
  "use strict";

  function init() {
    var input = document.querySelector(".search-box input[type='search']");
    var results = document.querySelector(".search-box .search-results");
    if (!input || !results) return;

    var base = input.getAttribute("data-base") || "";
    if (base && !/\/$/.test(base)) base += "/";

    var index = null;
    var activeIdx = -1;

    function loadIndex() {
      if (index !== null) return Promise.resolve(index);
      return fetch(base + "assets/search-index.json", { cache: "no-store" })
        .then(function (r) { return r.json(); })
        .then(function (data) {
          index = data.map(function (entry) {
            var hay = (entry.title + " " + (entry.section || "") + " " + (entry.kw || "")).toLowerCase();
            return Object.assign({}, entry, { _hay: hay });
          });
          return index;
        })
        .catch(function () {
          results.innerHTML = "<div class='no-match'>Search index not available.</div>";
          results.classList.add("open");
          return [];
        });
    }

    function render(query) {
      activeIdx = -1;
      var q = (query || "").trim().toLowerCase();
      if (!q) {
        results.innerHTML = "";
        results.classList.remove("open");
        return;
      }
      var terms = q.split(/\s+/);
      var matches = index.filter(function (e) {
        return terms.every(function (t) { return e._hay.indexOf(t) !== -1; });
      }).slice(0, 8);

      if (!matches.length) {
        results.innerHTML = "<div class='no-match'>No matches for &ldquo;" + escapeHtml(query) + "&rdquo;.</div>";
        results.classList.add("open");
        return;
      }

      results.innerHTML = matches.map(function (e) {
        return (
          "<a href='" + base + e.url + "'>" +
          "<div class='result-title'>" + escapeHtml(e.title) + "</div>" +
          "<div class='result-section'>" + escapeHtml(e.section || "") + "</div>" +
          "</a>"
        );
      }).join("");
      results.classList.add("open");
    }

    function escapeHtml(s) {
      return String(s).replace(/[&<>"']/g, function (c) {
        return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
      });
    }

    function setActive(delta) {
      var items = results.querySelectorAll("a");
      if (!items.length) return;
      activeIdx = (activeIdx + delta + items.length) % items.length;
      items.forEach(function (el, i) { el.classList.toggle("active", i === activeIdx); });
      var el = items[activeIdx];
      var box = results.getBoundingClientRect();
      var ebox = el.getBoundingClientRect();
      if (ebox.bottom > box.bottom) results.scrollTop += ebox.bottom - box.bottom;
      else if (ebox.top < box.top) results.scrollTop -= box.top - ebox.top;
    }

    input.addEventListener("focus", function () {
      loadIndex().then(function () { render(input.value); });
    });

    input.addEventListener("input", function () {
      loadIndex().then(function () { render(input.value); });
    });

    input.addEventListener("keydown", function (ev) {
      if (ev.key === "ArrowDown") { ev.preventDefault(); setActive(+1); }
      else if (ev.key === "ArrowUp") { ev.preventDefault(); setActive(-1); }
      else if (ev.key === "Enter") {
        var items = results.querySelectorAll("a");
        if (activeIdx >= 0 && items[activeIdx]) {
          ev.preventDefault();
          window.location.href = items[activeIdx].getAttribute("href");
        } else if (items.length) {
          ev.preventDefault();
          window.location.href = items[0].getAttribute("href");
        }
      } else if (ev.key === "Escape") {
        results.classList.remove("open");
        input.blur();
      }
    });

    document.addEventListener("click", function (ev) {
      if (!ev.target.closest(".search-box")) {
        results.classList.remove("open");
      }
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
