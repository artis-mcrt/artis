// Fill the version selector of the sidebar with the versions that Read the Docs holds.
// The reader then goes from the documentation of the main branch to the documentation of a release.
// https://docs.readthedocs.com/platform/stable/addons.html describes the event and the data.
function fillVersionSelector(versions) {
  const selector = document.querySelector(".version-selector");
  if (selector === null) {
    return;
  }
  const list = selector.querySelector("select");
  // a second call must give the same list, and not the list twice
  list.replaceChildren();
  // a hidden version and a build of a pull request stay out of the active versions, so the list needs the current one
  const shown = versions.active.some((version) => version.slug === versions.current.slug)
    ? versions.active
    : [versions.current, ...versions.active];
  for (const version of shown) {
    const entry = document.createElement("option");
    entry.value = version.urls.documentation;
    entry.text = version.slug;
    entry.selected = version.slug === versions.current.slug;
    list.add(entry);
  }
  list.onchange = () => {
    window.location.href = list.value;
  };
  selector.hidden = false;
}

document.addEventListener("readthedocs-addons-data-ready", (event) => {
  fillVersionSelector(event.detail.data().versions);
});

// Read the Docs loads its script with "async", so the event can come before this file runs. The data of the event
// stays in window.ReadTheDocsEventData, which gives the versions to a listener that starts too late.
if (window.ReadTheDocsEventData !== undefined) {
  fillVersionSelector(window.ReadTheDocsEventData.data().versions);
}
