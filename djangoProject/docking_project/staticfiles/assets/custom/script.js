


document.addEventListener('DOMContentLoaded', (event) => {
    console.log("hello");

    // Get the popup
    var popup = document.getElementById("popup3");

    // Get the button that opens the popup
    var btn = document.getElementById("viewDetailsBtn3");

    // Get the <span> element that closes the popup
    var closePopup = document.getElementById("closePopup3");

    // When the user clicks the button, open the popup 
    btn.onclick = function() {
        popup.style.display = "block";
    }

    // When the user clicks on <span> (x), close the popup
    closePopup.onclick = function() {
        popup.style.display = "none";
    }

    // When the user clicks anywhere outside of the popup, close it
    window.onclick = function(event) {
        if (event.target == popup) {
            popup.style.display = "none";
        }
    }
});

function generateLigand() {
    const virusName = document.getElementById('virusName').value;
    const virusCharacteristics = document.getElementById('virusCharacteristics').value;

    // Implement your ligand generation logic here
    console.log('Virus Name:', virusName);
    console.log('Virus Characteristics:', virusCharacteristics);

    // Example of what you might do next:
    alert(`Generating ligand for ${virusName} with characteristics: ${virusCharacteristics}`);
}
const nav = document.querySelector(".nav");
const navVideo = document.querySelector(".nav-menu-bg-vid video");

window.onscroll = () => {
    let scrollY = window.pageYOffset;
    document.querySelectorAll(".nav-dark").forEach((section) => {
        const sectionTop = section.offsetTop - 50;
        const sectionHeight = section.offsetHeight;

    });
    let navBlurOffset = 400;
    if (window.innerWidth < 490) {
        navBlurOffset = 80;
    }
};
const initTargetCards = function() {
    if (window.innerWidth > 991) {
        const cards = document.querySelectorAll(".science-targets-item");
        document.querySelector(".science-targets-list").style.width =
            (cards[0].clientWidth + 24) * cards.length -
            document.querySelector(".container").clientWidth +
            "px";
    }
}
initTargetCards();
window.addEventListener("resize", initTargetCards);
! function(o, c) {
    var n = c.documentElement,
        t = " w-mod-";
    n.className += t + "js", ("ontouchstart" in o || o.DocumentTouch && c instanceof DocumentTouch) && (n.className += t + "touch")
}(window, document);